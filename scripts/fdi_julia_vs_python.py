#!/usr/bin/env python

"""
"""

# 64 cells, 10 time steps
results = {
	'python': {
		'loop': 2*60 + 12.104,
		'vect': 0.899,
	},
	'julia': {
		'loop': 7*60 + 34.832,
		'vect': 1*60 + 46.731,
	},

}

import matplotlib
import numpy 

from matplotlib import pylab
font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
width = 0.6

count = 0
colors = ['orange', 'cyan', 'red', 'blue']
legends = []
for lang in ('python', 'julia'):
	for meth in ('loop', 'vect'):
		speed = results['python']['loop']/results[lang][meth]
		legends.append(lang + ' ' + meth)
		r = ax.bar([count*width], [speed], width, color=colors[count], log=True)
		count += 1

ax.set_ylabel('Speedup')
ax.set_xlabel('')
ax.set_title('Julia vs Python')
#ax.set_xticks(inds)
#ax.set_xticklabels(['{}'.format(n) for n in nth])
#ax.legend(list(platforms) + ['GPU'], loc=4)
ax.legend(legends)
plt.show()

