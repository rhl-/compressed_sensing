import cvxpy as cp
import mosek
import mosek.fusion
from   mosek.fusion import *
import numpy as np
import cvxopt
import common
import cmath
import sys
import scipy
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
	   return Matrix.sparse(m,p,i,j,v)

	# Make mosek environment
	def RPSD(A,y,n):
		with Model("RPSD_SOLVE") as M:
			X = M.variable("X", Domain.inPSDCone( n))
			obj = Expr.sum(X.diag())
			M.objective(ObjectiveSense.Minimize, obj)

			(i,j) = np.nonzero(A)
			v = A[i,j]
			(m,p) = A.shape
			B = Matrix.sparse(m,p,i,j,v)

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

			(i,j) = np.nonzero(A)
			v = A[i,j]
			(m,p) = A.shape
			B = Matrix.sparse(m,p,i,j,v)

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
			bigMat = M.variable(Domain.inPSDCone(2*n))

			U = bigMat.slice([0, 0], [n, n])
			V = bigMat.slice([n, 0], [2*n, n])

			obj = Expr.sum(U.diag())

			M.objective(ObjectiveSense.Minimize, obj)

			M.constraint(Expr.add(V, Variable.transpose(V)),Domain.equalsTo(0.0))

			B1 = get_matrix(A1)
			B2 = get_matrix(A2)
			# (i,j) = np.nonzero(A1)
			# v = A1[i,j]
			# (m,p) = A1.shape
			# B1 = Matrix.sparse(m,p,i,j,v)

			# (i,j) = np.nonzero(A2)
			# v = A2[i,j]
			# (m,p) = A2.shape
			# B2 = Matrix.sparse(m,p,i,j,v)

			RealPart = Expr.sub(Expr.mul(B1,Variable.reshape(U,n*n)),Expr.mul(B2,Variable.reshape(V,n*n)))
			ImagPart = Expr.add(Expr.mul(B1,Variable.reshape(V,n*n)),Expr.mul(B2,Variable.reshape(U,n*n)))

			M.constraint(RealPart,Domain.equalsTo(y1))
			M.constraint(ImagPart,Domain.equalsTo(y2))
			# M.writeTask("temp.opf")

			M.solve()
			try:


				U1 = np.matrix(U.level()).reshape(n,n)
				V1 = np.matrix(V.level()).reshape(n,n)
				# A = np.vstack((np.hstack((U1,-V1)),np.hstack((V1,U1))))

				# scipy.io.savemat('junk.mat', mdict={'A':A})

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

			bigMatPlus = M.variable(Domain.inPSDCone(2*n))
			bigMatMinus = M.variable(Domain.inPSDCone(2*n))


			Xp_r = bigMatPlus.slice([0, 0], [n, n])
			Xp_i = bigMatPlus.slice([n, 0], [2*n, n])

			Xm_r = bigMatMinus.slice([0, 0], [n, n])
			Xm_i = bigMatMinus.slice([n, 0], [2*n, n])

			M.constraint(Expr.add(Xp_i, Variable.transpose(Xp_i)),Domain.equalsTo(0.0))
			M.constraint(Expr.add(Xm_i, Variable.transpose(Xm_i)),Domain.equalsTo(0.0))


			obj = Expr.add(Expr.sum(Xp_r.diag()), Expr.sum(Xm_r.diag()))
			M.objective(ObjectiveSense.Minimize, obj)



			B1 = get_matrix(A1)
			B2 = get_matrix(A2)
			# (i,j) = np.nonzero(A1)
			# v = A1[i,j]
			# (m,p) = A1.shape
			# B1 = Matrix.sparse(m,p,i,j,v)

			# (i,j) = np.nonzero(A2)
			# v = A2[i,j]
			# (m,p) = A2.shape
			# B2 = Matrix.sparse(m,p,i,j,v)

			RealPart_p = Expr.sub(Expr.mul(B1,Variable.reshape(Xp_r,n*n)),Expr.mul(B2,Variable.reshape(Xp_i,n*n)))
			RealPart_m = Expr.sub(Expr.mul(B1,Variable.reshape(Xm_r,n*n)),Expr.mul(B2,Variable.reshape(Xm_i,n*n)))

			RealPart = Expr.sub(RealPart_p, RealPart_m)

			ImagPart_p = Expr.add(Expr.mul(B1,Variable.reshape(Xp_i,n*n)),Expr.mul(B2,Variable.reshape(Xp_r,n*n)))
			ImagPart_m = Expr.add(Expr.mul(B1,Variable.reshape(Xm_i,n*n)),Expr.mul(B2,Variable.reshape(Xm_r,n*n)))

			ImagPart = Expr.sub(ImagPart_p, ImagPart_m)

			M.constraint(RealPart,Domain.equalsTo(y1))
			M.constraint(ImagPart,Domain.equalsTo(y2))

			M.solve()

			try:
				Xp_numer = np.matrix(Xp_r.level()).reshape(n,n) + 1j * np.matrix(Xp_i.level()).reshape(n,n)
				Xm_numer = np.matrix(Xm_r.level()).reshape(n,n) + 1j * np.matrix(Xm_i.level()).reshape(n,n)

				x = Xp_numer-Xm_numer
				obj = np.sum(Xp_r.diag().level()) + np.sum(Xm_r.diag().level())

				#print np.linalg.norm(A*x.reshape(n*n,1) - y.reshape(m,1))
				return (obj, x,  cp.OPTIMAL)
			except:
				print "BAD"
				return (0, 0, 666)


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
