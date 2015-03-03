import mc_driver
from timer import Timer
t = Timer()

for i in xrange(3,8):
	print "k = %s"%(i)
	t.tic()
	mc_driver.run_trials(k=i,nMC=3,num_processes=1)
	t.toc()
