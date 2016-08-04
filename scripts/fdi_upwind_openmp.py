#!/usr/bin/env python

import matplotlib
import numpy 

from matplotlib import pylab
font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

# num cells = 128
# num steps = 100
# 1 thread

results = {
    'idefix': {
        'cpu': 'Intel i7',
        'clock GHz': 2.2,
        'compiler': 'clang-7.3',
        'num cores': 4,
        'time s': { 
            1: 35.6,
            2: float('nan'),
            4: float('nan'),
            8: float('nan'),
            16: float('nan'),
        },
    },
	'abraracourcix': {
		'cpu': 'Intel i7',
        'clock GHz': 2.2,
		'compiler': 'gnu-6.1',
		'num cores': 4,
		'time s': {	
			1: 27.7,
			2: 16.8,
			4: 9.5,
			8: 7.7,
			16: 7.7,
		},
	},
	'niwa-1007520': {
        'cpu': 'Intel Xeon',
        'clock GHz': 3.7,
        'compiler': 'gnu-4.8.5',
        'num cores': 4, # to check
        'time s': {
            1: 18.7,
            2: 12.7,
            4: 8.26,
            8: 5.46,
            16: 4.90,
            32: 4.56,
            64: 4.35,
        },
    },
    'niwa-1007520 pgi': {
        'cpu': 'Intel Xeon',
        'clock GHz': 3.7,
        'compiler': 'pgi16-3',
        'num cores': 4, # to check
        'time s': {
            1: 23.203,
            2: 15.149,
            4: 10.159,
            8: 6.619,
            16: 7.752,
            32: 12.52,
            64: float('nan'),
        },
    },
    'fitzroy': {
        'cpu': 'Power6',
        'clock GHz': 4.7,
        'compiler': 'xlcpp-12.1.0.13',
        'num cores': 32,
        'time s': {
            1: 141.7,
            2: 106.4,
            4: float('nan'),
            8: float('nan'),
            16: float('nan'),
            32: float('nan'),
            64: float('nan'),
        },
    },
    'pan-sb': {
        'cpu': 'Intel Xeon',
        'clock GHz': 3.7,
        'compiler': 'gnu-4.9.2',
        'num cores': 8, # to check
        'time s': {
            1: 25.1,
            2: 14.5,
            4: 9.8,
            8: 6.9,
            16: 4.9,
            32: float('nan'),
            64: float('nan'),
        },
    },
    'niwa-1007838': {
        'cpu': 'Intel i5',
        'clock GHz': 2.3,
        'compiler': 'gnu-4.8.5',
        'num cores': 2, 
        'time s': {
            1: 33.9,
            2: 24.1,
            4: 23.1,
            8: float('nan'),
        },
    },

}

# plot speed vs clock (1 thread)
clock = []
speed = []
cpu = []
for name, data in results.items():
    cpu.append(data['cpu'])
    clock.append(data['clock GHz'])
    speed.append(1./float(data['time s'][1]))

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.scatter(clock, speed, color='r', s=100)
for i in range(len(clock)):
    ax.annotate(cpu[i], (clock[i], speed[i]))
plt.xlim((0, 5))
plt.ylim((0, 0.08))
plt.xlabel('Clock GHz')
plt.ylabel('Speed 1/s')
plt.title('Speed vs clock frequency (128 cells 100 steps)')
#plt.show()

# how size affects the speed
ncells = [64, 128, 256, 512]
times_s = [2.32, 18.9, 150.2, 1191]

fig, ax = plt.subplots()
width = 0.2

count = 0
colors = ['c', 'b', 'orange']
#platforms = ('abraracourcix', 'pan-sb', 'niwa-1007520', 'niwa-1007520 pgi')
platforms = ('niwa-1007520', 'niwa-1007520 pgi')
for c in platforms:
    nth = sorted(results[c]['time s'])
    inds = [i + count*width for i in range(len(nth))]
    speed = [results['niwa-1007520']['time s'][1]/results[c]['time s'][n] for n in nth]
    r = ax.bar(inds, speed, width, color=colors[count])
    count += 1

# gpu speed
gpu_time_s = 4.85
ns = [1]
ids = [i + count*width for i in range(len(ns))]
speed = [results['niwa-1007520']['time s'][1]/gpu_time_s]
r = ax.bar(ids, speed, width, color='r')
ax.set_ylabel('Speedup')
ax.set_xlabel('number of CPU threads')
ax.set_title('OpenACC vs OpenMP')
ax.set_xticks(inds)
ax.set_xticklabels(['{}'.format(n) for n in nth])
#ax.legend(list(platforms) + ['GPU'], loc=4)
ax.legend(['OpenMP gcc', 'OpenMP pgi', 'OpenACC pgi'], loc=4)
plt.show()

