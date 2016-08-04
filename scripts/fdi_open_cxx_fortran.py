#!/usr/bin/env python

import matplotlib
import numpy 

from matplotlib import pylab
font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

import matplotlib.pyplot as plt

# timings in seconds
# niwa-1007520
results = {
	'gcc': {
        'c++': {
          1: 18.722,
          2: 12.441,
          4: 8.896,
          8: 5.592,
          16: 4.830,
        },
        'fortran': {
          1: 14.303,
          2: 6.882,
          4: 3.917,
          8: 3.064,
          16: 3.305,
        },
    },
    'pgi': {
        'c++': {
          1: 22.697,
          2: 14.948,
          4: 9.881,
          8: 6.490,
          16: 7.956,
        },
        'fortran': {
          1: 14.551,
          2: 9.458,
          4: 5.721,
          8: 5.005,
          16: 5.580,
        },
    },
}

fig, ax = plt.subplots()
width = 0.2

count = 0
colors = ['indigo', 'firebrick', 'c', 'orange']
platforms = ('gcc', 'pgi')
compilers = ('c++', 'fortran')
legs = []
for p in platforms:
    for c in compilers:
        legs.append(p + ' ' + c)
        nth = results[p][c].keys()
        nth.sort()
        inds = numpy.array(range(len(nth))) + count*width
        speed = [results['gcc']['c++'][1]/results[p][c][n] for n in nth]
        r = ax.bar(inds, speed, width, color=colors[count])
        count += 1

ax.set_ylabel('speedup')
ax.set_xlabel('number of OpenMP threads')
ax.set_title('C++ vs Fortran')
ax.set_xticks(inds)
ax.set_xticklabels(['{}'.format(n) for n in nth])
ax.legend(legs, loc=2)
plt.show()

