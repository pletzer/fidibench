#!/usr/bin/env python

import matplotlib
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
    'pan-sb': {
        'cpu': 'Intel Xeon',
        'clock GHz': 3.7,
        'compiler': 'gnu-4.8.2',
        'num cores': 8, # to check
        'time s': {
            1: float('nan'),
            2: float('nan'),
            4: float('nan'),
            8: float('nan'),
            16: float('nan'),
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
plt.xlim((0, 4))
plt.ylim((0, 0.08))
plt.xlabel('Clock GHz')
plt.ylabel('Speed 1/s')
plt.title('Speed vs clock frequency (128 cells 100 steps)')
plt.show()