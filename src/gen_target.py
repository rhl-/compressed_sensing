# Name:    gen_target.py
# Authors: Ryan Lewis & Victor Minden
# Purpose: A Python module to generate target matrices to recover via nuclear norm minimization.

import common
import numpy as np
import cmath

def getTargetSample(n, r, target):
	# getTargetSample(n, r, target)
	#
	# Input:
	# n       - generated matrix is of size n by n
	# r       - rank of matrix is r <= n
	# target  - an enum type taking on one of the following (enum) values:
    #           RPSD, RSYM, HPSD, HERM
	#
	# Output:
	# X0   -  an n by n array of rank r of type target

	def rand_diag():
		return np.diag(np.sign(np.random.normal(size=r)))
	#qr
	# LinAlgError :
		#If factoring fails.
	def makeRPSD():
		q, R = np.linalg.qr( np.random.normal(size=(n,r)),mode='reduced')
		#TODO if( q.shape[ 1] is not r)
			#ERROR
		return np.dot(q,q.transpose())

	def makeRSYM():
		q, R = np.linalg.qr( np.random.normal(size=(n,r)),mode='reduced')
		D = np.dot(rand_diag(),q.transpose())
		return np.dot(q, D)

	def makeHPSD():
		q = np.random.normal(size=(n,r)) + 1j*np.random.normal(size=(n,r))
		q, R = np.linalg.qr( q, mode='reduced')
		return np.dot(q, q.conj().transpose())

	def makeHERM():
		q = np.random.normal(size=(n,r)) + 1j*np.random.normal(size=(n,r))
		q, R = np.linalg.qr( q, mode='reduced')
		D = np.dot(rand_diag(),q.conj().transpose())
		return np.dot(q, D)

	if target == common.TARGET_TYPES.RPSD:
		return makeRPSD()
	elif target == common.TARGET_TYPES.RSYM:
		return makeRSYM()
	elif target == common.TARGET_TYPES.HPSD:
		return makeHPSD()
	elif target == common.TARGET_TYPES.HERM:
		return makeHERM()
	else:
		raise TypeError('Target type string does not match any in enum!')

