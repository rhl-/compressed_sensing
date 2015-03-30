#!/usr/bin/python
import numpy as np
import matplotlib as mpl
import sys
import os
import common

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)
ensembles = [common.ENSEMBLE_NAMES[0], common.ENSEMBLE_NAMES[2], common.ENSEMBLE_NAMES[8]]
targets = [common.TARGET_NAMES[2], common.TARGET_NAMES[3]]
#targets = [common.TARGET_NAMES[0], common.target_names[ 3]]
mpl.use('agg')
import matplotlib.pyplot as plt

for filename in sys.argv[1:]:
	outputfilestem = filename.split('.')[0]
	# outputfilename = os.path.splitext(filename)
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


	bp = ax.boxplot(time_plot_data.values(),patch_artist=True)
	ax.set_xticklabels(time_plot_data.keys())
	ax.set_ylabel('Time (s)')
	ax.set_title('Timing')

		## change outline color, fill color and linewidth of the boxes
	for box in bp['boxes']:
	    # change outline color
	    box.set( color='#7570b3', linewidth=2)
	    # change fill color
	    box.set( facecolor = '#1b9e77' )

	## change color and linewidth of the whiskers
	for whisker in bp['whiskers']:
	    whisker.set(color='#7570b3', linewidth=2)

	## change color and linewidth of the caps
	for cap in bp['caps']:
	    cap.set(color='#7570b3', linewidth=2)

	## change color and linewidth of the medians
	for median in bp['medians']:
	    median.set(color='#b2df8a', linewidth=2)

	## change the style of fliers and their fill
	for flier in bp['fliers']:
	    flier.set(marker='o', color='#e7298a', alpha=0.5)

	plt.show()
	fig.savefig(outputfilestem + '_time.pdf', bbox_inches='tight')


	fig = plt.figure(2, figsize=(20,20))
	ax = fig.add_subplot(111)

	bp = ax.boxplot(mem_plot_data.values(),patch_artist=True)
	ax.set_xticklabels(mem_plot_data.keys())
	ax.set_ylabel('TB')
	ax.set_title('Memory usage')

		## change outline color, fill color and linewidth of the boxes
	for box in bp['boxes']:
	    # change outline color
	    box.set( color='#7570b3', linewidth=2)
	    # change fill color
	    box.set( facecolor = '#1b9e77' )

	## change color and linewidth of the whiskers
	for whisker in bp['whiskers']:
	    whisker.set(color='#7570b3', linewidth=2)

	## change color and linewidth of the caps
	for cap in bp['caps']:
	    cap.set(color='#7570b3', linewidth=2)

	## change color and linewidth of the medians
	for median in bp['medians']:
	    median.set(color='#b2df8a', linewidth=2)

	## change the style of fliers and their fill
	for flier in bp['fliers']:
	    flier.set(marker='o', color='#e7298a', alpha=0.5)

	plt.show()
	fig.savefig(outputfilestem + '_mem.pdf', bbox_inches='tight')
