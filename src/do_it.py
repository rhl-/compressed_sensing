import mc_driver
import common
from timer import Timer
t = Timer()

#ensemble, target, nmc, RPSD, RDIRAC
#mc_driver.run_rank_trial(m=1024,n=64,ranks=xrange(5,16),ensemble=0,target=0,nMC=25)
mc_driver.run_trials(ms=range(100,4001,200),k=6,nMC=20)
#print "k = 5"
#t.tic()
#mc_driver.run_trials([(32**2)/8],k=5)
#t.toc()
#
#
#print "k = 6"
#t.tic()
#mc_driver.run_trials([(64**2)/8],k=6)
#t.toc()
#
#
#print "k = 5"
#t.tic()
#mc_driver.run_trials([(32**2)/4],k=5)
#t.toc()
#
#
#print "k = 6"
#t.tic()
#mc_driver.run_trials([(64**2)/4],k=6)
#t.toc()







# print "k = 6"
# t.tic()
# mc_driver.run_trials(k=6)
# t.toc()

# print "k = 7"
# t.tic()
# mc_driver.run_trials(k=7)
# t.toc()
