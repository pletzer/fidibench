from matplotlib import rcParams, pylab

rcParams['font.size'] = 14

"""
ml cray-libsci_acc/18.06.1 craype-accel-nvidia60 PrgEnv-cray/1.0.4 cuda92/blas/9.2.88 cuda92/toolkit/9.2.88
time srun -n 1 -t 20 -p gpu --gres=gpu:1 ./upwindOMP5Cxx -numSteps 100 -numCells 256
for n in 36 36 18 18 10 10 4 4 2 2 1; do echo $n; export OMP_NUM_THREADS=$n; export OMP_PROC_BIND=true export OMP_PLACES=cores; time srun --nodes=1 --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS ./upwind2Cxx -numSteps 100 -numCells 256; done
"""


data = {
	'mahuika-cray': {
		'nthreads': [36, 18, 9, 4, 2, 1],
		'seconds' : [4.6, 7.9, 11.6, 24.5, 46.3, 60.7],
		'color': 'r',
		'style': '--'
		},


}

d = data['mahuika-cray']
pylab.plot(d['nthreads'], d['seconds'], d['color'] + d['style'])
pylab.plot([min(d['nthreads']), max(d['nthreads'])], [7.6, 7.6], 'c-')
pylab.legend(['OpenMP-CPU', 'OpenMP-GPU'])
pylab.xlabel('num CPU threads')
pylab.ylabel('time [s]')
pylab.title('fidibench upwind C++ 100 steps/256^3 cells')
pylab.axis([1, 40, 0, 80])
pylab.show()



