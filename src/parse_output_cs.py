#!/usr/bin/python
import numpy as np
import matplotlib as mpl
import sys
import common

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)
ensembles = [common.ENSEMBLE_NAMES[0], common.ENSEMBLE_NAMES[2], common.ENSEMBLE_NAMES[8]]
targets = [common.TARGET_NAMES[0], common.TARGET_NAMES[1]]
#targets = [common.TARGET_NAMES[0], common.target_names[ 3]]
mpl.use('agg')
import matplotlib.pyplot as plt

for filename in sys.argv[1:]:
	f = open(filename)
	lines = f.readlines()
	time_plot_data={ (x,y):[] for x in targets for y in ensembles }
	mem_plot_data={ (x,y):[] for x in targets for y in ensembles }
	n = ""
	for line_str in lines:
		line = line_str.split(',')
		if( line[1] in targets and line[2] in ensembles):
			key = (line[1], line[2])
			time_plot_data[key].append( float(line[-2]))
			mem_plot_data[key].append( float(line[-1]))
			n = line[3]	
	fig = plt.figure(1, figsize=(20,20))
	ax = fig.add_subplot(111)
	bp = ax.boxplot(time_plot_data.values())
	ax.set_xticklabels(time_plot_data.keys())
	plt.show()
	fig.savefig('figtest.pdf', bbox_inches='tight')
