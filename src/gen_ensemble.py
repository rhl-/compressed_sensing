# Name:    gen_ensemble.py
# Authors: Ryan Lewis & Victor Minden
# Purpose: A Python module to generate sampling ensembles
#	  for matrix recovery via nuclear norm minimization.

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
        #           ENTRY,  PERM,   RSPERM, CSPERM, RGPERM,
	# 	    CGPERM, RDIRAC, CDIRAC, RGAUSS, CGAUSS
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

	# A random permutation
	def make_PERM():
		A = np.zeros((m,n**2))
		for i in xrange(m):
			permi = np.random.permutation(n)
			Ai = np.zeros((n,n))
			for (j,k) in enumerate(permi):
				Ai[j,k] = 1.0
			A[i,:] = Ai.flatten()
		return A/np.sqrt(n)

	# Add a sign to the permutation
	def make_RSPERM():
		A = make_PERM()
		# THIS RANDOM CALL IS DENSE AND EXPENSIVE FOR BIG N
		A[np.where(np.random.rand(m,n**2) > 0.5)] *= -1.0
		return A

	# Modulate with a complex phase
	def make_CSPERM():
		B = make_RSPERM()
		A = np.zeros((m,n**2),dtype=complex) + B
		A[np.where(np.random.rand(m,n**2) > 0.5)] *= 1.0j
		return A

	# Modulate with a Gaussian
	def make_RGPERM():
		A = make_PERM()
		return np.multiply(A, np.random.normal(size=(m,n**2)))

	def make_CGPERM():
		#Initially A is unscaled by 1/sqrt(n)
		A = make_PERM()
		return np.multiply((1.0/np.sqrt(2.0))*A,(np.random.normal(size=(m,n**2))) + 1j*np.random.normal(size=(m,n**2)))

	def generate_kronecker_sequence():
		p = np.random.permutation(n**2)
		if (target == common.TARGET_TYPES.RPSD) or (target == common.TARGET_TYPES.RSYM):
			p = [ i for i in p if (np.base_repr( i, 4).count('2')%2 ==0)]
		return p[:m]

	def dirac_matrix( basis):
		assert len(basis) == 4
		A = np.zeros((m,n**2))
		J = np.log2(n)
		assert J == round(J)
		for (i,val) in enumerate(generate_kronecker_sequence()):
			seq = [int(x) for x in list( np.base_repr(val, 4))]
			seq = [0]*(J-len(seq)) + seq
			Ai = np.matrix('1.0')
			for j in seq:
				Ai = np.kron(Ai,basis[j])
			A[i,:] = Ai.flatten()
		return A/np.sqrt(n)

	def make_RDIRAC():
		I = np.matrix('1.0 0.0; 0.0 1.0')
		wx = np.matrix('0.0 1.0; 1.0 0.0')
		wy = np.matrix('0.0 -1.0; 1.0 0.0')
		wz = np.matrix('1.0 0.0; 0.0 -1.0')
        	basis =(I,wx,wy,wz)
        	return dirac_matrix(basis);

	def make_CDIRAC():
		I = np.matrix('1.0 0.0; 0.0 1.0')
		sx = np.matrix('0.0 1.0; 1.0 0.0')
		sy = np.matrix('0.0 -sqrt(-1); sqrt(-1) 0.0')
		sz = np.matrix('1 0.0; 0.0 -1.0')
        	basis =(I,sx,sy,sz)
        	return dirac_matrix(basis);

	def make_RGAUSS():
		return (1.0/n)*np.random.normal(size=(m,n*n))

	def make_CGAUSS():
        	return (1.0/(np.sqrt(2)*n)) * (np.random.normal(size=(m,n**2)) + 1j*np.random.normal(size=(m,n**2)))

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
		return make_CGAUSS()
	else:
		raise TypeError('Ensemble type string does not match any in enum!')
