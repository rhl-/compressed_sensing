# Name:    gen_ensemble.py
# Authors: Ryan Lewis & Victor Minden
# Purpose: A Python module to generate sampling ensembles for matrix recovery via nuclear norm minimization.

import common

def matrix_sample(m, n, meas, mat):
	# Input:
	# m    - number of measurements
	# n    - dimension of matrixes we will be "sensing" with the ensemble
	# meas - an enum type taking on one of the following (enum) values:
    #        ENTRY, PERM, RSPERM, CSPERM, RGPERM, CGPERM, RDIRAC, CDIRAC, RGAUSS, CGAUSS
    # mat  - an enum type taking on one of the following (enum) values:
    #        RPSD, RSYM, HPSD, HERM
    #
	# Output:
	# A    - array of size m by n^2, containing the sampling matrices as row vectors


	def make_ENTRY():
		return 0
	def make_PERM():
		return 1
	def make_RSPERM():
		return 2
	def make_CSPERM():
		return 3
	def make_RGPERM():
		return 4
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


