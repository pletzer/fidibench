#!/usr/bin/env python

from matplotlib import pylab

# num cells = 128
# num steps = 100
# 1 thread

results {
	'abraracourcix': {
		'cpu': 'Intel i7 2.2GHz',
		'compiler': 'gnu-6.1',
		'num cores': 4,
		'time_s': {	
			1: 27.7,
			2: 16.8,
			4: 9.5,
			8: 7.7,
			16: 7.7,
		}
	},
	'idefix': {},
	'pan-wm': {},
	'pan-sb': {},
	'niwa-1007520': {
                'cpu': 'Intel Xeon 3.7GHz',
                'compiler': 'gnu-4.8.5',
                'num cores': 4, # to check
                'time_s': {
                        1: 18.7,
                        2: 12.7,
                        4: 8.26,
                        8: 5.46,
                        16: 4.90,
                        32: 4.56,
                        64: 4.35,
        },
        'niwa-1007838': {
                'cpu': 'Intel i5 2.3GHz',
                'compiler': 'gnu-4.8.5',
                'num cores': 2, 
                'time_s': {
                        1: ,
                        2: ,
                        4: ,
                        8: ,
        }

	
}
