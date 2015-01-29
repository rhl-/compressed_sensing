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

	def makeRPSD():
		pass
	def makeRSYM():
		pass
	def makeHPSD():
		pass
	def makeHERM():
		pass

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

