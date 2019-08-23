from matplotlib import rcParams, pylab

rcParams['font.size'] = 14

"""
export OMP_NUM_THREADS=2; export OMP_PROC_BIND=true; export OMP_PLACES=cores; srun --nodes=1 --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --mem=5g map --profile ./upwindCxx -numSteps 10 -numCells 512

mahuika
for n in 36 18 9 4 2 1; do echo $n; export OMP_NUM_THREADS=$n; export OMP_PROC_BIND=true; export OMP_PLACES=cores; srun --nodes=1 --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --mem=5g time ./upwindCxx -numSteps 10 -numCells 512; done

maui
for n in 40 20 10 4 2 1; do echo $n; export OMP_NUM_THREADS=$n; export OMP_PROC_BIND=true; export OMP_PLACES=cores; srun --nodes=1 --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --mem=5g time ./upwindCxx -numSteps 10 -numCells 512; done
"""

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
	pylab.axis([1, 40, 0, 150])
	pylab.show()

data = {
	'mahuika-cray': {
		'nthreads': [36, 18, 9, 4, 2, 1],
		'seconds-old' : [131, 102, 102, 86, 130, 189],
		'seconds-new' : [6.7, 10.1, 21.3, 36.3, 74.4, 145],
		'color': 'r',
		'style': '--'
		},
	'mahuika-intel': {
		'nthreads': [36, 18, 9, 4, 2, 1],
		'seconds-old' : [31, 68, 78, 109, 303, 182],
		'seconds-new' : [6.2, 9.8, 17.5, 35.8, 70, 139],
		'color': 'g',
		'style': '--'
		},
	'mahuika-gnu': {
		'nthreads': [36, 18, 9, 4, 2, 1],
		'seconds-old' : [ 22, 63, 84, 132, 197, 177,], # [45, 180, 236, 361, 447, 807],
		'seconds-new' : [ 6.5, 10.2, 20.6, 37.5, 73.0, 134], # [18, 32, 58, 124, 242, 483],
		'color': 'b',
		'style': '--'
		},

	'maui-cray': {
		'nthreads': [40, 20, 10, 4, 2, 1],
		'seconds-old' : [7.3, 12, 18, 39, 71, 124],
		'seconds-new' : [6.5, 9.9, 13.9, 27.8, 51.0, 98.4],
		'color': 'r',
		'style': '-'
		},

	'maui-intel': {
		'nthreads': [40, 20, 10, 4, 2, 1],
		'seconds-old' : [27, 33, 81, 99, 183, 126],
		'seconds-new' : [5.8, 8.5, 13.0, 27.5, 49.5, 94.7],
		'color': 'g',
		'style': '-'
		},

	'maui-gnu': {
		'nthreads': [40, 20, 10, 4, 2, 1],
		'seconds-old' : [32, 32, 26, 85, 137, 115],
		'seconds-new' : [6.6, 9.1, 13.4, 27.1, 47.8, 90.8],
		'color': 'b',
		'style': '-'
		},

}

plotNew(keys=['mahuika-' + p for p in ('cray', 'intel', 'gnu')])
plotNew(keys=['maui-' + p for p in ('cray', 'intel', 'gnu')])
plotNew(keys=['mahuika-' + p for p in ('cray', 'intel', 'gnu')] + ['maui-' + p for p in ('cray', 'intel', 'gnu')])


