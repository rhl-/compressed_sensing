# Name:    gen_ensemble.py
# Authors: Ryan Lewis & Victor Minden
# Purpose: A Python module to generate sampling ensembles for matrix recovery via nuclear norm minimization.

import common
import numpy as np
import cmath

def getEnsembleSample(m, n, meas, target):
	# getEnsembleSample(m, n, meas, target)
	#
	# Input:
	# m       - number of measurements
	# n       - dimension of matrixes we will be "sensing" with the ensemble
	# meas 	  - an enum type taking on one of the following (enum) values:
    #           ENTRY, PERM, RSPERM, CSPERM, RGPERM, CGPERM, RDIRAC, CDIRAC, RGAUSS, CGAUSS
    # target  - an enum type taking on one of the following (enum) values:
    #           RPSD, RSYM, HPSD, HERM
    #
	# Output:
	# A       - array of size m by n^2, containing the sampling matrices as row vectors


	def make_ENTRY():
		# population is the linear indices of the strict upper triangle of an n by n matrix
		population = np.flatnonzero(np.triu(np.ones((n,n))))
		# idxs is a set of random indices without repetition
		idxs = np.random.choice(population,size=m,replace=False)
		A = np.zeros((m,n**2))
		for (i,j) in enumerate(idxs):
			A[i,j] = 1.0
		return A

	def make_PERM():
		# Each matrix is a random permutation scaled by 1/sqrt(n)
		A = np.zeros((m,n**2))
		for i in xrange(m):
			permi = np.random.permutation(n)
			Ai = np.zeros((n,n))
			for (j,k) in enumerate(permi):
				Ai[j,k] = 1.0 / np.sqrt(n)
			A[i,:] = Ai.flatten()
		return A

	def make_RSPERM():
		# Add a sign to the permutation
		A = make_PERM()
		A[np.where(np.random.rand(m,n**2) > 0.5)] *= -1.0
		return A

	def make_CSPERM():
		# Modulate with a complex phase
		B = make_RSPERM()
		A = np.zeros((m,n**2),dtype=complex)
		A += B
		A[np.where(np.random.rand(m,n**2) > 0.5)] *= 1.0j
		return A

	def make_RGPERM():
		# Modulate with a Gaussian
		A = make_RSPERM()
		A = np.multiply(A, np.random.normal(size=(m,n**2)))
		return A

	def make_CGPERM():
		return 5
	def make_RDIRAC():
		return 6
	def make_CDIRAC():
		return 7
	def make_RGAUSS():
		return 8
	def make_CGAUSS():
		return 9


	# Each subfunction will define the returned matrix A,
	# so just call the appropriate subfunction and then return A

	if meas == common.ENSEMBLE_TYPES.ENTRY:
		return make_ENTRY()
	elif meas == common.ENSEMBLE_TYPES.PERM:
		return make_PERM()
	elif meas == common.ENSEMBLE_TYPES.RSPERM:
		return make_RSPERM()
	elif meas == common.ENSEMBLE_TYPES.CSPERM:
		return make_CSPERM()
	elif meas == common.ENSEMBLE_TYPES.RGPERM:
		return make_RGPERM()
	elif meas == common.ENSEMBLE_TYPES.CGPERM:
		return make_CGPERM()
	elif meas == common.ENSEMBLE_TYPES.RDIRAC:
		return make_RDIRAC()
	elif meas == common.ENSEMBLE_TYPES.CDIRAC:
		return make_CDIRAC()
	elif meas == common.ENSEMBLE_TYPES.RGAUSS:
		return make_RGAUSS()
	elif meas == common.ENSEMBLE_TYPES.CGAUSS:
		return make_GAUSS()
	else:
		raise TypeError('Ensemble type string does not match any in enum!')


