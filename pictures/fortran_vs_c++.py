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
intel_fortran_times = [6.96, 9.13, 13.46, 22.38, 40.88, 1*60+18.1]
intel_cxx_times = [6.5, 10.5, 18.7, 34.9, 1*60+7.8, 2*60+13.5]


pylab.figure(1)
pylab.plot(nthreads_mahuika, intel_cxx_times, 'r-')
pylab.plot(nthreads_mahuika, intel_fortran_times, 'b-')
pylab.legend(['C++ intel 19.1.0.166', 'Fortran intel 19.1.0.166'])
pylab.xlabel('number of OpenMP threads')
pylab.ylabel('execution time [s]')
pylab.title('mahuika numCells: 512 numSteps: 10')
pylab.plot(nthreads_mahuika, intel_cxx_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_mahuika, intel_fortran_times, 'ko', markerfacecolor='None')

pylab.show()
