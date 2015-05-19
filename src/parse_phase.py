#!/usr/bin/python
import numpy as np
import matplotlib as mpl
import sys
import os
import common

#mpl.use('agg')
import matplotlib.pyplot as plt

filename = sys.argv[1]
ens_idx  = int(sys.argv[2])
targ_idx = int(sys.argv[3])

ens = common.ENSEMBLE_NAMES[ens_idx]
targ = common.TARGET_NAMES[targ_idx]


f = open(filename)
lines = f.readlines()

err0 = []
rank = []
meas = []

for line_str in lines:
	line = line_str.split(',')
	if ( line[2] == targ and line[3] == ens):
		err0.append( float(line[8]))
		rank.append( int(line[5]))
		meas.append( int(line[6]))


points = {}

ranks = []
mcs   = []
for (r,m) in zip(rank,meas):
	points[(r,m)] = [0,0]

for (e,r,m) in zip(err0,rank,meas):
	temp = points[(r,m)]
	temp[0] += e
	temp[1] += 1

ranks = []
meass = []
colors = []

for key in points:
	ranks.append(float(key[0]) / 64)
	meass.append(float(key[1]) / (64 ** 2))
	colors.append(float(points[key][0]) / points[key][1])



# for key in points
# for key in points:
# 	ranks.append(key)
# 	mcs.append(float(points[key][0]) / points[key][1])

print ranks
print meass
print colors

plt.scatter(ranks, meass, s=500,c=colors)
# plt.xticks(np.arange(1, max(ranks), 2))
plt.xlabel('r/N')
plt.ylabel('m/N^2')
plt.title('N=64, Ensemble=%s, Target=%s'%(common.ENSEMBLE_NAMES[ens_idx], common.TARGET_NAMES[targ_idx]))
plt.show()


		#time_plot_data[key].append( float(line[-2]))
		# mem_plot_data[key].append( float(line[-1]))
		#n = line[3]
# fig = plt.figure(1, figsize=(20,20))
# ax = fig.add_subplot(111)


	# bp = ax.boxplot(time_plot_data.values(),patch_artist=True)
	# ax.set_xticklabels(time_plot_data.keys())
	# ax.set_ylabel('Time (s)')
	# ax.set_title('Timing')

	# 	## change outline color, fill color and linewidth of the boxes
	# for box in bp['boxes']:
	#     # change outline color
	#     box.set( color='#7570b3', linewidth=2)
	#     # change fill color
	#     box.set( facecolor = '#1b9e77' )

	# ## change color and linewidth of the whiskers
	# for whisker in bp['whiskers']:
	#     whisker.set(color='#7570b3', linewidth=2)

	# ## change color and linewidth of the caps
	# for cap in bp['caps']:
	#     cap.set(color='#7570b3', linewidth=2)

	# ## change color and linewidth of the medians
	# for median in bp['medians']:
	#     median.set(color='#b2df8a', linewidth=2)

	# ## change the style of fliers and their fill
	# for flier in bp['fliers']:
	#     flier.set(marker='o', color='#e7298a', alpha=0.5)

	# plt.show()
	# fig.savefig(outputfilestem + '_time.pdf', bbox_inches='tight')


	# fig = plt.figure(2, figsize=(20,20))
	# ax = fig.add_subplot(111)

	# bp = ax.boxplot(mem_plot_data.values(),patch_artist=True)
	# ax.set_xticklabels(mem_plot_data.keys())
	# ax.set_ylabel('TB')
	# ax.set_title('Memory usage')

	# 	## change outline color, fill color and linewidth of the boxes
	# for box in bp['boxes']:
	#     # change outline color
	#     box.set( color='#7570b3', linewidth=2)
	#     # change fill color
	#     box.set( facecolor = '#1b9e77' )

	# ## change color and linewidth of the whiskers
	# for whisker in bp['whiskers']:
	#     whisker.set(color='#7570b3', linewidth=2)

	# ## change color and linewidth of the caps
	# for cap in bp['caps']:
	#     cap.set(color='#7570b3', linewidth=2)

	# ## change color and linewidth of the medians
	# for median in bp['medians']:
	#     median.set(color='#b2df8a', linewidth=2)

	# ## change the style of fliers and their fill
	# for flier in bp['fliers']:
	#     flier.set(marker='o', color='#e7298a', alpha=0.5)

	# plt.show()
	# fig.savefig(outputfilestem + '_mem.pdf', bbox_inches='tight')
