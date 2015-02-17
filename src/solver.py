import cvxpy as cp
import mosek
import mosek.fusion
from   mosek.fusion import *
import numpy as np
import cvxopt
import common
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
				return (np.sum(X.diag().level()), x,  cp.OPTIMAL)
			except:
				return (0, 0, 666)


	# def RPSD_cp(A,y,n):
	# 	X = cp.Semidef(n)
	# 	objective = cp.Minimize( cp.norm( X, "nuc"))
	# 	constraints = [A*cp.vec(X)==y]
	# 	# Construct and solve the cvx_py problem
	# 	prob = cp.Problem( objective, constraints)
	# 	prob.solve(verbose=False)

	# 	return (prob.value, X.value, prob.status)

	# def RPSD(A,y,n):
	# 	(val1,X1,stat1) = RPSD_mos(A,y,n)
	# 	(val2,X2,stat2) = RPSD_cp(A,y,n)
	# 	print np.linalg.norm(X1-X2,'fro')
	# 	print val1
	# 	print val2

	def RSYM(A,y,n):
		X = cp.Variable(n,n)
		objective = cp.Minimize( cp.norm( X, "nuc"))
		constraints = [A*cp.vec(X)==y, X == X.T]
		prob =  cp.Problem( objective, constraints)
		prob.solve(verbose=False)
		return (prob.value, X.value, prob.status)

	def HPSD(A,y,n):
		with Model("HPSD_SOLVE") as M:
			# W = [Real(X), -Imag(X); Imag(X), Real(X)]
			U = M.variable(NDSet(n,n), Domain.unbounded())
			V = M.variable(NDSet(n,n), Domain.unbounded())
			obj = Expr.sum(U.diag())

			M.objective(ObjectiveSense.Minimize, obj)

			row1 = Expr.hstack(U,Expr.neg(V));
			row2 = Expr.hstack(V,U);


			M.constraint(Expr.vstack(row1,row2), Domain.inPSDCone(2*n))

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
				U1 = np.array(U.level()).reshape(n,n)
				V1 = np.array(V.level()).reshape(n,n)
				print np.linalg.norm(V1)

				x = U1+cmath.sqrt(-1) * V1

				return (2*np.sum(U.diag().level()), x,  cp.OPTIMAL)
			except:
				return (0, 0, 666)



	def HERM(A,y,n):
		U = cp.variable(n,n)
		V = cp.variable(n,n)

		A1 = np.real(A)
		A2 = np.imag(A)
		y1 = np.real(y)
		y2 = np.imag(y)

		objective = cp.Minimize( cp.norm( U, "nuc"))
		constraints = [A1*cp.vec(U)==y1, A2*cp.vec(V) == y2, U == U.T, V == -V.T]
		prob = cp.Problem( objective, constraints)
		prob.solve(verbose=False)
		return (prob.value, U.value + sqrt(-1) * V.value, prob.status)

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
