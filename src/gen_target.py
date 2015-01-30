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
		return np.diag(np.sign(np.random.normal(size=(n,r))))
	#qr
	# LinAlgError :
		#If factoring fails.
	def makeRPSD():
		q = np.linalg.qr( np.random.normal(size=(n,r)),mode='economic')
		#TODO if( q.shape[ 1] is not r)
			#ERROR
		return q*q.transpose()

	def makeRSYM():
		q = np.linalg.qr( np.random.normal(size=(n,r)),mode='economic')
		return q*rand_diag()*q.transpose() 

	def makeHPSD():
		q = np.random.normal(size=(n,r)) + 1j*np.random.normal(size=(n,r))
		q = np.linalg.qr( q, mode='economic')
		return q*q.transpose()

	def makeHERM():
		q = np.random.normal(size=(n,r)) + 1j*np.random.normal(size=(n,r))
		q = np.linalg.qr( q, mode='economic')
		return q*rand_diag()*q.transpose() 

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

