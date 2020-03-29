import matplotlib
import numpy as np
import matplotlib.pyplot as plt
font = {'family' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
from matplotlib import pylab

# for nth in 36 24 18 12 4 2 1; do export OMP_NUM_THREADS=$nth; time ./upwind/cxx/upwindCxx -numCells 512 -numSteps 10 ; done

# upwindCxx -numCells 512 -numSteps 10
nthreads_mahuika = [32, 16, 8, 4, 2, 1]

# Release
intel_times = [6.5, 10.5, 18.7, 34.9, 1*60+7.8, 2*60+13.5]
gnu_times = [39.6, 59.8, 1*60+16.3, 1*60+26.4, 3*60+11.6, 2*60+52.8]
gnuO3_times = [6.6, 10.4, 18.4, 34.4, 1*60+6.4, 2*60+10.7]
cray_times = [6.8, 11.0, 19.2, 37.1, 1*60+9.9, 2*60+17.5]


pylab.figure(1)
pylab.plot(nthreads_mahuika, intel_times, 'r-')
pylab.plot(nthreads_mahuika, gnu_times, 'g--')
pylab.plot(nthreads_mahuika, gnuO3_times, 'g--')
pylab.plot(nthreads_mahuika, cray_times, 'b-')
pylab.legend(['intel 19.1.0.166', 'gnu 9.2.0', 'gnu -O3 9.2.0', 'cray 8.7.7'])
pylab.xlabel('number of OpenMP threads')
pylab.ylabel('execution time [s]')
pylab.title('mahuika upwindCxx -numCells 512 -numSteps 10')
pylab.plot(nthreads_mahuika, intel_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_mahuika, gnu_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_mahuika, gnuO3_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_mahuika, cray_times, 'ko', markerfacecolor='None')

pylab.show()
