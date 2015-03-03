import cvxpy as cp
import mosek
import mosek.fusion
from   mosek.fusion import *
import numpy as np
import cvxopt
import common
import cmath
import sys
# CVX algorithm for nuclear norm minimization on square matrices
# Input:
# A : the measurement matrix A. Preconditions # rows of A = n^2
# y : the measurements
# n : first shape dimension of X
# type : 'RPSD','RSYM','HPSD','HERM'
#
# Output:
# X: the CVX solution.
# val: NaN if the algorithm diverged
def solve(A,y,n,t):
	def get_matrix( A):
	   (m,p) = A.shape
	   nnz = np.count_nonzero(A)
	   if( nnz == m*p):
	   	return DenseMatrix( A)
	   
	   (i,j) = np.nonzero(A)
	   v = A[i,j]
	   (m,p) = A.shape
	   B = Matrix.sparse(m,p,i,j,v)
			
	# Make mosek environment
	def RPSD(A,y,n):
		with Model("RPSD_SOLVE") as M:
			X = M.variable("X", Domain.inPSDCone( n))
			obj = Expr.sum(X.diag())
			M.objective(ObjectiveSense.Minimize, obj)

			B = get_matrix( A)
	
			M.constraint(Expr.mul(B,Variable.reshape(X,n*n)),Domain.equalsTo(y))

			M.solve()
			try:
				x = np.array(X.level()).reshape(n,n)
				#print np.linalg.norm(x-x.T)
				return (np.sum(X.diag().level()), x,  cp.OPTIMAL)
			except:
				return (0, 0, 666)


	def RSYM(A,y,n):
		with Model("RSYM_SOLVE") as M:
			Xp = M.variable("Xp", Domain.inPSDCone( n))
			Xm = M.variable("Xm", Domain.inPSDCone( n))

			# Add the two traces
			obj = Expr.add(Expr.sum(Xp.diag()), Expr.sum(Xm.diag()))
			M.objective(ObjectiveSense.Minimize, obj)

			B = get_matrix( A)
			
			#X = Xp - Xm
			X = Expr.sub(Xp,Xm)
			M.constraint(Expr.mul(B,Expr.reshape(X,n*n)),Domain.equalsTo(y))

			M.solve()
			try:
				x = np.array(Xp.level()).reshape(n,n) - np.array(Xm.level()).reshape(n,n)
				return (np.sum(Xp.diag().level()) + np.sum(Xm.diag().level()), x,  cp.OPTIMAL)
			except:
				print "BAD"
				return (0, 0, 666)


	def HPSD(A,y,n):
		A1 = np.real(A)
		A2 = np.imag(A)
		y1 = np.real(y)
		y2 = np.imag(y)

		with Model("HPSD_SOLVE") as M:
			U = M.variable("U", Domain.inPSDCone(n))
			V = M.variable(NDSet(n,n), Domain.unbounded())

			obj = Expr.sum(U.diag())

			M.objective(ObjectiveSense.Minimize, obj)

			row1 = Expr.hstack(U,Expr.mul(-1.0,V));
			row2 = Expr.hstack(V,U);

			M.constraint(Expr.vstack(row1,row2), Domain.inPSDCone(2*n))

			M.constraint(Expr.add(V, Variable.transpose(V)),Domain.equalsTo(0.0))
		
			(i,j) = np.nonzero(A1)
			v = A1[i,j]
			(m,p) = A1.shape
			B1 = Matrix.sparse(m,p,i,j,v)

			(i,j) = np.nonzero(A2)
			v = A2[i,j]
			(m,p) = A2.shape
			B2 = Matrix.sparse(m,p,i,j,v)

			RealPart = Expr.sub(Expr.mul(B1,Variable.reshape(U,n*n)),Expr.mul(B2,Variable.reshape(V,n*n)))
			ImagPart = Expr.add(Expr.mul(B1,Variable.reshape(V,n*n)),Expr.mul(B2,Variable.reshape(U,n*n)))

			M.constraint(RealPart,Domain.equalsTo(y1))
			M.constraint(ImagPart,Domain.equalsTo(y2))
			# M.writeTask("temp.opf")

			M.solve()
			try:
				U1 = np.matrix(U.level()).reshape(n,n)
				V1 = np.matrix(V.level()).reshape(n,n)
				#A = np.vstack((np.hstack((U1,-V1)),np.hstack((V1,U1))))
				#print A

				x = U1 + 1j* V1
				obj = 2*np.sum(U.diag().level())

				# print A*x.reshape(n*n,1)
				# print y
				#print np.linalg.norm(A*x.reshape(n*n,1) - y.reshape(m,1))
				return (obj, x,  cp.OPTIMAL)
			except:

				print "BAD"
				return (0, 0, 666)



	def HERM(A,y,n):
		A1 = np.real(A)
		A2 = np.imag(A)
		y1 = np.real(y)
		y2 = np.imag(y)
		with Model("HERM_SOLVE") as M:
			# W = [Real(X), -Imag(X); Imag(X), Real(X)]
			# Real part is in PSD cone, imaginary part is skew-symmetric
			Xp_r = M.variable("Xp_r", Domain.inPSDCone( n))
			Xp_i = M.variable(NDSet(n,n), Domain.unbounded())

			Xm_r = M.variable("Xm_r", Domain.inPSDCone( n))
			Xm_i = M.variable(NDSet(n,n), Domain.unbounded())

			#Skew symmetric on imaginary parts
			# M.constraint(Expr.add(Xp_i, Variable.transpose(Xp_i)),Domain.equalsTo(0.0))
			# M.constraint(Expr.add(Xm_i, Variable.transpose(Xm_i)),Domain.equalsTo(0.0))

			# Trace (Real(Xp) + Real(Xm))
			obj = Expr.add(Expr.sum(Xp_r.diag()), Expr.sum(Xm_r.diag()))
			M.objective(ObjectiveSense.Minimize, obj)

			ext_Xp_row1 = Expr.hstack(Xp_r,Expr.mul(-1.0,Xp_i))
			ext_Xp_row2 = Expr.hstack(Xp_i,Xp_r)

			ext = Expr.vstack(ext_Xp_row1,ext_Xp_row2)
			#print ext
			M.constraint(ext, Domain.inPSDCone(2*n))


			ext_Xm_row1 = Expr.hstack(Xm_r,Expr.mul(-1.0,Xm_i))
			ext_Xm_row2 = Expr.hstack(Xm_i,Xm_r)

			M.constraint(Expr.vstack(ext_Xm_row1,ext_Xm_row2), Domain.inPSDCone(2*n))

			(i,j) = np.nonzero(A1)
			v = A1[i,j]
			(m,p) = A1.shape
			B1 = Matrix.sparse(m,p,i,j,v)

			(i,j) = np.nonzero(A2)
			v = A2[i,j]
			(m,p) = A2.shape
			B2 = Matrix.sparse(m,p,i,j,v)

			RealPart_p = Expr.sub(Expr.mul(B1,Variable.reshape(Xp_r,n*n)),Expr.mul(B2,Variable.reshape(Xp_i,n*n)))
			RealPart_m = Expr.sub(Expr.mul(B1,Variable.reshape(Xm_r,n*n)),Expr.mul(B2,Variable.reshape(Xm_i,n*n)))

			RealPart = Expr.sub(RealPart_p, RealPart_m)

			ImagPart_p = Expr.add(Expr.mul(B1,Variable.reshape(Xp_i,n*n)),Expr.mul(B2,Variable.reshape(Xp_r,n*n)))
			ImagPart_m = Expr.add(Expr.mul(B1,Variable.reshape(Xm_i,n*n)),Expr.mul(B2,Variable.reshape(Xm_r,n*n)))

			ImagPart = Expr.sub(ImagPart_p, ImagPart_m)

			M.constraint(RealPart,Domain.equalsTo(y1))
			M.constraint(ImagPart,Domain.equalsTo(y2))
			#M.setLogHandler(sys.stdout)
			M.solve()
			# M.writeTask("temp.opf")
			# print U.level()
			try:
				Xp_numer = np.matrix(Xp_r.level()).reshape(n,n) + 1j * np.matrix(Xp_i.level()).reshape(n,n)
				Xm_numer = np.matrix(Xm_r.level()).reshape(n,n) + 1j * np.matrix(Xm_i.level()).reshape(n,n)
				#A = np.vstack((np.hstack((U1,-V1)),np.hstack((V1,U1))))
				#print A

				x = Xp_numer-Xm_numer
				obj = np.sum(Xp_r.diag().level()) + np.sum(Xm_r.diag().level())
				#t = np.matrix(Xp_r.level()).reshape(n,n)
				#U,s,V = np.linalg.svd(t)
				#print s
				#res = A*x.reshape(n*n,1) - y.reshape(m,1)
				# print res.shape
				# print np.linalg.norm(res)
				# print A*x.reshape(n*n,1)
				# print y
				#print np.linalg.norm(A*x.reshape(n*n,1) - y.reshape(m,1))
				return (obj, x,  cp.OPTIMAL)
			except:
				print "BAD"
				return (0, 0, 666)
		# A1 = np.real(A)
		# A2 = np.imag(A)
		# y1 = np.real(y)
		# y2 = np.imag(y)
		# with Model("HERM_SOLVE") as M:
		# 	# W = [Real(X), -Imag(X); Imag(X), Real(X)]
		# 	U = M.variable(NDSet(n,n), Domain.unbounded())
		# 	V = M.variable(NDSet(n,n), Domain.unbounded())

		# 	Z1 = M.variable("Z1", Domain.inPSDCone( n))
		# 	Z2 = M.variable(NDSet(n,n), Domain.unbounded())

		# 	obj = Expr.sum(Z1.diag())

		# 	M.objective(ObjectiveSense.Minimize, obj)

		# 	row1 = Expr.hstack(U,Expr.mul(-1.0,V));
		# 	row2 = Expr.hstack(V,U);

		# 	B1 = Expr.vstack(Expr.hstack(Z1,U), Expr.hstack(U,Z1))
		# 	B2 = Expr.vstack(Expr.hstack(Z2,V), Expr.hstack(V,Z2))


		# 	BLOCK = Expr.vstack( Expr.hstack(B1, Expr.neg(B2)), Expr.hstack(B2, B1))

		# 	M.constraint(BLOCK, Domain.inPSDCone(4*n))


		# 	M.constraint(Expr.add(V, Variable.transpose(V)),Domain.equalsTo(0.0))
		# 	M.constraint(Expr.sub(U, Variable.transpose(U)),Domain.equalsTo(0.0))
		# 	M.constraint(Expr.sub(Z1, Variable.transpose(Z1)),Domain.equalsTo(0.0))
		# 	M.constraint(Expr.add(Z2, Variable.transpose(Z2)),Domain.equalsTo(0.0))

		# 	(i,j) = np.nonzero(A1)
		# 	v = A1[i,j]
		# 	(m,p) = A1.shape
		# 	B1 = Matrix.sparse(m,p,i,j,v)

		# 	(i,j) = np.nonzero(A2)
		# 	v = A2[i,j]
		# 	(m,p) = A2.shape
		# 	B2 = Matrix.sparse(m,p,i,j,v)

		# 	RealPart = Expr.sub(Expr.mul(B1,Variable.reshape(U,n*n)),Expr.mul(B2,Variable.reshape(V,n*n)))
		# 	ImagPart = Expr.add(Expr.mul(B1,Variable.reshape(V,n*n)),Expr.mul(B2,Variable.reshape(U,n*n)))

		# 	M.constraint(RealPart,Domain.equalsTo(y1))
		# 	M.constraint(ImagPart,Domain.equalsTo(y2))
		# 	# M.setLogHandler(sys.stdout)
		# 	M.solve()
		# 	# M.writeTask("temp.opf")
		# 	# print U.level()
		# 	try:
		# 		U1 = np.matrix(U.level()).reshape(n,n)
		# 		V1 = np.matrix(V.level()).reshape(n,n)
		# 		#A = np.vstack((np.hstack((U1,-V1)),np.hstack((V1,U1))))
		# 		#print A

		# 		x = U1 + 1j* V1
		# 		obj = np.sum(Z1.diag().level())

		# 		# print A*x.reshape(n*n,1)
		# 		# print y
		# 		#print np.linalg.norm(A*x.reshape(n*n,1) - y.reshape(m,1))
		# 		return (obj, x,  cp.OPTIMAL)
		# 	except:

		# 		print "BAD"
		# 		return (0, 0, 666)
	# 	U = cp.variable(n,n)
	# 	V = cp.variable(n,n)

	# 	A1 = np.real(A)
	# 	A2 = np.imag(A)
	# 	y1 = np.real(y)
	# 	y2 = np.imag(y)

	# 	objective = cp.Minimize( cp.norm( U, "nuc"))
	# 	constraints = [A1*cp.vec(U)==y1, A2*cp.vec(V) == y2, U == U.T, V == -V.T]
	# 	prob = cp.Problem( objective, constraints)
	# 	prob.solve(verbose=False)
	# 	return (prob.value, U.value + sqrt(-1) * V.value, prob.status)

	if t == common.TARGET_TYPES.RPSD:
		return RPSD(A,y,n)
	elif t == common.TARGET_TYPES.RSYM:
		return RSYM(A,y,n)
	elif t == common.TARGET_TYPES.HPSD:
		return HPSD(A,y,n)
	elif t == common.TARGET_TYPES.HERM:
		return HERM(A,y,n)
	else:
		raise
