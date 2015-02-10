import cvxpy as cp
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
	def RPSD(A,y,n):

		X = cp.Semidef(n)

		objective = cp.Minimize( cp.norm( X, "nuc"))

		# Construct and solve the cvx_py problem
		prob = cp.Problem( objective, [A*cp.vec(X)==y])
		prob.solve(verbose=True)

		return (prob.value, X.value)

	def RSYM(A,y,n):

		X = cp.Variable(n,n)
		objective = cp.Minimize( cp.norm( X, "nuc"))
		prob =  cp.Problem( objective, [A*cp.vec(X)==y, X == X.T])
		prob.solve(verbose=True)
		return (prob.value, X.value)

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
		prob = cp.Problem( objective, [A1*cp.vec(U)==y1, A2*cp.vec(V) == y2, W == Semidef(n)])
		prob.solve(verbose=True)
		return (prob.value, U.value + sqrt(-1) * V.value)

	def HERM(A,y,n):
		U = cp.variable(n,n)
		V = cp.variable(n,n)

		A1 = np.real(A)
		A2 = np.imag(A)
		y1 = np.real(y)
		y2 = np.imag(y)

		objective = cp.Minimize( cp.norm( U, "nuc"))
		prob = cp.Problem( objective, [A1*cp.vec(U)==y1, A2*cp.vec(V) == y2, U == U.T, V == -V.T])
		prob.solve(verbose=True)
		return (prob.value, U.value + sqrt(-1) * V.value)

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
