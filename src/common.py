# Name:    common.py
# Authors: Ryan Lewis & Victor Minden
# Purpose: Contains utility routines for project on matrix recovery via nuclear norm minimization

def enum(*sequential, **named):
	# An enum type with example usage as:
	#	>>> Numbers = enum('ZERO', 'ONE', 'TWO')
	# 	>>> Numbers.ZERO
	# 	0
	# 	>>> Numbers.ONE
	# 	1
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

ENSEMBLE_NAMES = ['ENTRY', 'PERM', 'RSPERM', 'CSPERM', 'RGPERM', 'CGPERM', 'RDIRAC', 'CDIRAC', 'RGAUSS', 'CGAUSS']
ENSEMBLE_TYPES = enum('ENTRY', 'PERM', 'RSPERM', 'CSPERM', 'RGPERM', 'CGPERM', 'RDIRAC', 'CDIRAC', 'RGAUSS', 'CGAUSS')
TARGET_TYPES   = enum('RPSD', 'RSYM', 'HPSD', 'HERM')
TARGET_NAMES   = ['RPSD', 'RSYM', 'HPSD', 'HERM']
