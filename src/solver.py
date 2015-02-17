import cvxpy as cp
import mosek
import mosek.fusion
from   mosek.fusion import *
import numpy as np
import cvxopt
import common
import cmath
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
				print np.linalg.norm(x-x.T)
				return (np.sum(X.diag().level()), x,  cp.OPTIMAL)
			except:
				return (0, 0, 666)


	def RSYM(A,y,n):
		with Model("RSYM_SOLVE") as M:
			X = M.variable(NDSet( n,n), Domain.unbounded())
			U = M.variable("U", Domain.inPSDCone( n))

			obj = Expr.sum(U.diag())
			M.objective(ObjectiveSense.Minimize, obj)

			(i,j) = np.nonzero(A)
			v = A[i,j]
			(m,p) = A.shape
			B = Matrix.sparse(m,p,i,j,v)

			M.constraint(Expr.mul(B,Variable.reshape(X,n*n)),Domain.equalsTo(y))

			BLOCK = Variable.vstack( Variable.hstack(U, X), Variable.hstack(X, U))
			M.constraint( BLOCK , Domain.inPSDCone( 2*n) )

			M.constraint(Expr.sub(X,Variable.transpose(X)), Domain.equalsTo(0.0))

			M.solve()
			try:
				x = np.array(X.level()).reshape(n,n)
				u = np.array(U.level()).reshape(n,n)
				return (2*np.sum(U.diag().level()), x,  cp.OPTIMAL)
			except:
				return (0, 0, 666)


	def HPSD(A,y,n):
		A1 = np.real(A)
		A2 = np.imag(A)
		y1 = np.real(y)
		y2 = np.imag(y)

		with Model("HPSD_SOLVE") as M:
			# W = [Real(X), -Imag(X); Imag(X), Real(X)]
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
			U = M.variable(NDSet(n,n), Domain.unbounded())
			V = M.variable(NDSet(n,n), Domain.unbounded())

			Z = M.variable("Z", Domain.inPSDCone( 2*n))

			obj = Expr.sum(U.diag())
			obj = Expr.sum(U.diag())

			M.objective(ObjectiveSense.Minimize, obj)

			row1 = Expr.hstack(U,Expr.mul(-1.0,V));
			row2 = Expr.hstack(V,U);


			X = Expr.vstack(row1,row2)
			BLOCK = Variable.vstack( Variable.hstack(Z, X), Variable.hstack(X, Z))

			M.constraint(BLOCK, Domain.inPSDCone(4*n))


			M.constraint(Expr.add(V, Variable.transpose(V)),Domain.equalsTo(0.0))
			M.constraint(Expr.sub(U, Variable.transpose(U)),Domain.equalsTo(0.0))

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
