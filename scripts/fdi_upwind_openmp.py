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
	''
	
}