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
from multiprocessing import Process
import time
import datetime

class MCDriver:
	def __init__(self, k=3, nMC=10, filename=None):
		if(filename == None):
			date_str=str(datetime.datetime.now()).replace(" ","")
			filename = "matrix_completion_k_%s_nMc_%s_%s"%(k,nMC,date_str)
		np.random.seed(12181990)
		self.n = 2**k
		self.nMC = nMC
		self.filename = filename
		# target, ensemble, n, r, m, err_1, err_2, time_target, time_ensemble, time_solve
		self.formatstring = "%s,%s,%d,%d,%d,%f,%f,%f,%f,%f\n"


	def run_trials(self):
		# to sweep over m and r and generate a whole bunch of trials of
		# each and call solver on those while logging output

		for target in range(len(common.TARGET_NAMES)):
			for ensemble in range(len(common.ENSEMBLE_NAMES)):

				real_target = (target == common.TARGET_TYPES.RPSD or target == common.TARGET_TYPES.RSYM)
				complex_measurement = (common.ENSEMBLE_NAMES[ensemble][0] == 'C')

				if not (real_target and complex_measurement):
					# Choose random ranks
					if self.n < 32:
						l = 1
					else:
						l = 6
					rs = range(1,6,1) + randi(l,self.n/4,5).tolist() + randi(self.n/4+1,self.n/2,5).tolist() + randi(self.n/2+1,self.n,5).tolist()
					#print rs

					# # First, calculate upper bound on number of measurements
					# is_entry = (ensemble == common.ENSEMBLE_TYPES.ENTRY )
					# is_real_dirac = (common.ENSEMBLE_NAMES[ensemble][1:] == "DIRAC" and real_target)
					# if is_entry or is_real_dirac:
					# 	UB = self.n ** 2 / 2 + self.n/2
					# else:
					# 	UB = self.n ** 2
					UB = self.n ** 2 / 2 + self.n / 2

					ms_foreach_r = [[] for _ in rs]
					for (ms,r) in zip(ms_foreach_r,rs):
						rho = float(r) / self.n
						# Compute lower bound on number of measurements in a given test
						LB = max(np.rint((rho - rho ** 2 / 2) * self.n ** 2),1)
						ms.extend(randi(LB, UB, 10).tolist())

					proc_list = []
					for (ms, r) in zip(ms_foreach_r,rs):
						for m in ms:
							for trial in xrange(self.nMC):
								# Here is where we will spawn a new process
								p = Process(target=self.gen_problem, args=(m,r,target,ensemble,))
								proc_list.append(p)
								p.start()
								#self.gen_problem(m,r,target,ensemble)


	def gen_problem(self, m, r, target, meas):
		#print "%d %d %d %d" % (m,r,target,meas)
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


