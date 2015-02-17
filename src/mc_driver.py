# Name:    mc_driver.py
# Authors: Ryan Lewis & Victor Minden
# Purpose: A Python module that contains the Monte Carlo driver for the main experiment.
import common
from matrix_recovery_problem import *
import numpy as np # Might not need this
from numpy.random import randint as randi
from multiprocessing import Pool
import time
import datetime
import sys

def compute_task(n,nMC,filename,m,r,target,ensemble):
	#TODO: This can instead go ahead and call out for this..
	instance = problem_instance(n,nMC,filename)
	instance.solve_problem(m,r,target,ensemble)

def run_trials(k=3,nMC=10,filename=None):
	tic = time.time()
	if(filename == None):
		date_str=datetime.datetime.now().strftime("%b-%d-%y-%I:%M%p")
		filename = "OUTPUT_n_%s_nMc_%s_%s"%(2**k,nMC,date_str)
	np.random.seed(12181990)
	n = 2**k
	# to sweep over m and r and generate a whole bunch of trials of
	# each and call solver on those while logging output
	#p = Pool()
	for target in [2]:#range(len(common.TARGET_NAMES)):
		for ensemble in range(len(common.ENSEMBLE_NAMES)):
			#40 iterations of this loop it seems.
			real_target = (target == common.TARGET_TYPES.RPSD or target == common.TARGET_TYPES.RSYM)
			complex_measurement = (common.ENSEMBLE_NAMES[ensemble][0] == 'C')
			if not (real_target and complex_measurement):
				# Choose random ranks
				if n < 32:
					l = 1
				else:
					l = 6
				rs = range(1,6,1) + randi(l,n/4,5).tolist() + randi(n/4+1,n/2,5).tolist() + randi(n/2+1,n,5).tolist()
				#print rs

				# # First, calculate upper bound on number of measurements
				# is_entry = (ensemble == common.ENSEMBLE_TYPES.ENTRY )
				# is_real_dirac = (common.ENSEMBLE_NAMES[ensemble][1:] == "DIRAC" and real_target)
				# if is_entry or is_real_dirac:
				# 	UB = n ** 2 / 2 + n/2
				# else:
				# 	UB = n ** 2
				UB = n ** 2 / 2 + n / 2
				ms_foreach_r = [[] for _ in rs]
				for (ms,r) in zip(ms_foreach_r,rs):
					rho = float(r) / n
					# Compute lower bound on number of measurements in a given test
					LB = max(np.rint((rho - rho ** 2 / 2) * n ** 2),1)
					ms.extend(randi(LB, UB, 10).tolist())

				for (ms, r) in zip(ms_foreach_r,rs):
					for m in ms:
						for trial in xrange(nMC):
							#p.apply_async(compute_task, args=(n,nMC,filename,m,r,target,ensemble,))
							compute_task(n,nMC,filename,m,r,target,ensemble)


	toc = time.time()-tic
	print "Queued all jobs \t %f"%(toc)
	sys.stdout.write('Waiting...')
	sys.stdout.flush()
	tic=time.time()
	# p.close()
 #    	p.join()
	toc = time.time()-tic
	print " done."
	print"Spent %f seconds waiting for processes"%(toc)
