from matplotlib import rcParams, pylab

rcParams['font.size'] = 14


def plotOld(keys):
	for k in keys:
		d = data[k]
		pylab.plot(d['nthreads'], d['seconds-old'], 
			d['color'] + d['style'])
	pylab.legend(keys)
	pylab.xlabel('num threads')
	pylab.ylabel('time [s]')
	pylab.title('fidibench upwind C++ 10 steps/512^3 cells')
	pylab.axis([1, 40, 0, 800])
	pylab.show()

def plotNew(keys):
	for k in keys:
		d = data[k]
		pylab.plot(d['nthreads'], d['seconds-new'], 
			d['color'] + d['style'])
	pylab.legend(keys)
	pylab.xlabel('num threads')
	pylab.ylabel('time [s]')
	pylab.title('fidibench upwind C++ 10 steps/512^3 cells')
	pylab.axis([1, 40, 0, 800])
	pylab.show()

data = {
	'mahuika-cray': {
		'nthreads': [36, 18, 9, 4, 2, 1],
		'seconds-old' : [131, 102, 102, 86, 130, 189],
		'seconds-new' : [6.2, 11, 19, 36, 72, 142],
		'color': 'r',
		'style': '--'
		},
	'mahuika-intel': {
		'nthreads': [36, 18, 9, 4, 2, 1],
		'seconds-old' : [31, 68, 78, 109, 303, 182],
		'seconds-new' : [6.2, 10, 17, 36, 75, 138],
		'color': 'g',
		'style': '--'
		},
	'mahuika-gnu': {
		'nthreads': [36, 18, 9, 4, 2, 1],
		'seconds-old' : [ 22, 63, 84, 132, 197, 177,], # [45, 180, 236, 361, 447, 807],
		'seconds-new' : [ 6.1, 9.5, 17.3, 35.2, 68.6, 134], # [18, 32, 58, 124, 242, 483],
		'color': 'b',
		'style': '--'
		},

	'maui-cray': {
		'nthreads': [40, 20, 10, 4, 2, 1],
		'seconds-old' : [7.3, 12, 18, 39, 71, 124],
		'seconds-new' : [6.5, 9.8, 14, 27, 50, 96],
		'color': 'r',
		'style': '-'
		},

	'maui-intel': {
		'nthreads': [40, 20, 10, 4, 2, 1],
		'seconds-old' : [27, 33, 81, 99, 183, 126],
		'seconds-new' : [5.9, 8.6, 13, 28, 50, 96],
		'color': 'g',
		'style': '-'
		},

	'maui-gnu': {
		'nthreads': [40, 20, 10, 4, 2, 1],
		'seconds-old' : [32, 32, 26, 85, 137, 115],
		'seconds-new' : [6.4, 9.0, 13.2, 26.8, 48, 90],
		'color': 'b',
		'style': '-'
		},

}

plotOld(keys=['mahuika-' + p for p in ('cray', 'intel', 'gnu')])
plotOld(keys=['maui-' + p for p in ('cray', 'intel', 'gnu')])


