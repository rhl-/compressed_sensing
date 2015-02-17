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
			# B = DenseMatrix(A)
			# print B
			# print X
			# print "N is %s" % (n * n)
			# print y

			# print Variable.reshape(X,n*n)

			M.constraint(Expr.mul(B,Variable.reshape(X,n*n)),Domain.equalsTo(y))

			M.solve()
			try:
				x = np.array(X.level()).reshape(n,n)
				return (np.sum(X.diag().level()), x,  cp.OPTIMAL)
			except:
				return (0, 0, 666)
#		X = cp.Semidef(n)
#		objective = cp.Minimize( cp.norm( X, "nuc"))
#		constraints = [A*cp.vec(X)==y]
#		# Construct and solve the cvx_py problem
#		prob = cp.Problem( objective, constraints)
#		prob.solve(verbose=False)
#
#		return (prob.value, X.value, prob.status)

	def RSYM(A,y,n):

		X = cp.Variable(n,n)
		objective = cp.Minimize( cp.norm( X, "nuc"))
		constraints = [A*cp.vec(X)==y, X == X.T]
		prob =  cp.Problem( objective, constraints)
		prob.solve(verbose=False)
		return (prob.value, X.value, prob.status)

	def HPSD(A,y,n):
		# W = [Real(X), -Imag(X); Imag(X), Real(X)]
		U = cp.variable(n,n)
		V = cp.variable(n,n)
		W = cp.vstack(cp.hstack(U,-V), cp.hstack(V,U))

		A1 = np.real(A)
		A2 = np.imag(A)
		y1 = np.real(y)
		y2 = np.imag(y)

		objective = cp.Minimize( cp.norm( U, "nuc"))
		constraints = [A1*cp.vec(U)==y1, A2*cp.vec(V) == y2, W == Semidef(n)]
		prob = cp.Problem( objective, constraints)
		prob.solve(verbose=False)
		return (prob.value, U.value + sqrt(-1) * V.value, prob.status)

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
