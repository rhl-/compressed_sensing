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
		X = cp.Semidef( n)
		objective = cp.Minimize( cp.norm( X, "nuc"))
		return cp.Problem( objective, [A*cp.vec(X)==y])
	def RSYM(A,y,n):
		X = cp.symmetric( n)
		objective = cp.Minimize( cp.norm( X, "nuc"))
		return cp.Problem( objective, [A*cp.vec(X)==y])
	def HPSD(A,y,n):
		X = cp.hermitian( n)
		objective = cp.Minimize( cp.norm( X, "nuc"))
		return cp.Problem( objective, [A*cp.vec(X)==y, X == cp.hermitian_semidefinite( n)])
	def HERM(A,y,n):
		X = cp.hermitian( n)
		objective = cp.Minimize( cp.norm( X, "nuc"))
		return cp.Problem( objective, [A*cp.vec(X)==y])
	functions = [RPSD,RSYM,HPSD,HERM]
	options = dict(zip(common.TARGET_TYPES,functions))
	return options[ t]( A, y, n).solve()
