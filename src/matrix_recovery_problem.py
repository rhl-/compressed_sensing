# Name:    mc_driver.py
# Authors: Ryan Lewis & Victor Minden
# Purpose: A Python module that contains the Monte Carlo driver for the main experiment.
import cvxpy as cp
import common
import gen_ensemble
import gen_target
import numpy as np # Might not need this
from numpy.random import randint as randi
import solver
from lockfile import LockFile
import time
import datetime

class problem_instance:
	def __init__(self, n, nMC, filename):
		#TODO: Random seeds are perhaps better than victors birthday.
		np.random.seed(12181990)
		self.n = n
		self.nMC = nMC
		self.filename = filename
		# target, ensemble, n, r, m, err_1, err_2, time_target, time_ensemble, time_solve
		self.formatstring = "%s,%s,%d,%d,%d,%f,%f,%f,%f,%f\n"

	def solve_problem(self, m, r, target, meas):
		# generate a problem of rank r with m measurements
		tic = time.time()
		Xtrue = gen_target.getTargetSample(self.n, r, target)
		time_target = time.time()-tic
		tic = time.time()
		A = gen_ensemble.getEnsembleSample(m, self.n, meas, target)
		y = np.dot(A,Xtrue.flatten())
		time_ensemble = time.time()-tic
		tic = time.time()
		(opt_val, x_opt, status) = self.call_solver(A, y, self.n, target)
		time_solve = time.time()-tic
		if( status == cp.OPTIMAL):
			self.log_output(target, meas,x_opt,Xtrue, r, m,
time_target,time_ensemble,time_solve)
		else:
			print "Issue with solver! not logging..."
	def call_solver(self, A, y, n, target):
		return solver.solve(A, y, n, target)

	def log_output(self, target, meas, x_opt, Xtrue, r, m, t1,t2,t3):
		# write output to a file
		# we should use a different file for each trial and then just do
		# a "reduce".  This way any individual task isn't too long
		try:
			lock = LockFile(self.filename)
		except:
			raise "Could not open file"

		tol = 1e-2
		err0 = 0
		err1 = np.linalg.norm(Xtrue - x_opt,'fro') / np.linalg.norm(Xtrue)
		if err1 > tol:
			err0 = 1
		lock.acquire()
		with open(self.filename, 'a') as my_file:
			try:
				my_file.write(self.formatstring % (common.TARGET_NAMES[target],
common.ENSEMBLE_NAMES[meas], self.n, r, m, err1, err0, t1,t2,t3))
			except:
				print "Could not open file!"
		lock.release()
