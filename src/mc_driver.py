# Name:    mc_driver.py
# Authors: Ryan Lewis & Victor Minden
# Purpose: A Python module that contains the Monte Carlo driver for the main experiment.

import common
import gen_ensemble
import gen_target
import numpy as np # Might not need this

class MCDriver:
	def __init__(self, n=32, target=common.TARGET_TYPES.RPSD, nMC=10):
	# __init__(self, n, target, nMC)
	#
	# Input:
	# n       - target matrices are of size n by n
	# target  - an enum type taking on one of the following (enum) values:
	#           RPSD, RSYM, HPSD, HERM
	# nMC     - the number of Monte Carlo trials to use for each data point


    	# Need to set the number of measurements and the ranks to match what we 
	# did in Matlab
		pass
	def run_trials(self):
		# to sweep over m and r and generate a whole bunch of trials of 
		# each and call solver on those while logging output
		pass
	def gen_problem(self, r=self.n/2, m):
		# generate a problem of rank r with m measurements
		pass
	def call_solver(self):
		# either call cvx opt or fork to a C++ solver or something
		pass
	def log_output(self):
		# write output to a file
		# we should use a different file for each trial and then just do 
		# a "reduce".  This way any individual task isn't too long
		pass
