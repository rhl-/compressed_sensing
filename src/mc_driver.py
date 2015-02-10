# Name:    mc_driver.py
# Authors: Ryan Lewis & Victor Minden
# Purpose: A Python module that contains the Monte Carlo driver for the main experiment.

import common
import gen_ensemble
import gen_target
import numpy as np # Might not need this
import solver
class MCDriver:
	def __init__(self, n=32, target=common.TARGET_TYPES.RPSD, nMC=1):
	# __init__(self, n, target, nMC)
	#
	# Input:
	# n       - target matrices are of size n by n
	# target  - an enum type taking on one of the following (enum) values:
	#           RPSD, RSYM, HPSD, HERM
	# nMC     - the number of Monte Carlo trials to use for each data point


    	# Need to set the number of measurements and the ranks to match what we
	# did in Matlab
		self.n = n
		self.target = target
		self.nMC = nMC
		#build ranks and number of measurements, but for now use junk
		self.m = n
		self.r = 5
		pass
	def run_trials(self):
		# to sweep over m and r and generate a whole bunch of trials of
		# each and call solver on those while logging output

		#hack_test
		self.gen_problem(self.n, self.r)

	def gen_problem(self, m, r):
		# generate a problem of rank r with m measurements
		Xtrue = gen_target.getTargetSample(self.n, self.r, self.target)
		#hack
		meas = common.ENSEMBLE_TYPES.ENTRY
		A = gen_ensemble.getEnsembleSample(self.m, self.n, meas,self.target)
		y = np.dot(A,Xtrue.flatten())
		self.call_solver(A, y, self.n, self.target)
	# (A * vec(X))_i = y_i
	def call_solver(self, A, y, n, target):
		(opt_val, x_opt) = solver.solve(A, y, n, target)
		#HACK
		print opt_val

	def log_output(self):
		# write output to a file
		# we should use a different file for each trial and then just do
		# a "reduce".  This way any individual task isn't too long
		pass
