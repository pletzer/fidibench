import matplotlib
import numpy as np
import matplotlib.pyplot as plt
font = {'family' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
from matplotlib import pylab

# for nth in 36 24 18 12 4 2 1; do export OMP_NUM_THREADS=$nth; time ./upwind/cxx/upwindCxx -numCells 512 -numSteps 10 ; done

# upwindCxx -numCells 512 -numSteps 10
nthreads_mahuika = [36, 20, 10, 4, 2, 1]

# Release
intel_times = [ 38.898, 51.181, 1*60+17.584, 2*60+42.927, 2*60+59.922, 3*60+6.235]
gnu_times = [ 39.081, 55.064, 1*60+26.947, 2*60+54.221, 2*60+59.998, 3*60+6.689]
cray_times = [  58.714, 54.368, 1*60+1.994, 1*60+28.097, 2*60+25.811, 3*60+4.515]

#
# kupe -numCells 800 -numSteps 10
#
nthreads_kupe = [40, 20 , 10, 4, 2, 1]
intel_kupe_times = [ 1*60+7.612, 1*60+42.784, 1*60+31.348, 3*60+33.351, 4*60+35.903, 8*60+0.621]
intel_mahuika_times = [ 2*60+22.864, 3*60+33.783,  5*60+6.747, 9*60+41.842, 12*60+48.018, 11*60+58.291]
gnu71_kupe_times = [ 1*60+29.877, 3*60+11.809, 4*60+4.363, 6*60+32.641, 4*60+11.765, 7*60+21.261]
gnu61_kupe_times = [ 1*60+32.590, 3*60+13.742, 3*60+55.533, 6*60+56.419, 4*60+11.935, 7*60+20.378]
cray_kupe_times = [  31.925, 51.051, 1*60+16.631, 2*60+40.138, 4*60+59.781, 8*60+27.299]


pylab.figure(1)
pylab.plot(nthreads_mahuika, intel_times, 'r-')
pylab.plot(nthreads_mahuika, gnu_times, 'g-')
pylab.plot(nthreads_mahuika, cray_times, 'b-')
pylab.legend(['intel', 'gnu', 'cray'])
pylab.xlabel('number of OpenMP threads')
pylab.ylabel('execution time [s]')
pylab.title('mahuika upwindCxx -numCells 512 -numSteps 10')
pylab.plot(nthreads_mahuika, intel_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_mahuika, gnu_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_mahuika, cray_times, 'ko', markerfacecolor='None')

pylab.show()

pylab.figure(2)
pylab.plot(nthreads_kupe, intel_kupe_times, 'r-')
pylab.plot(nthreads_mahuika, intel_mahuika_times, 'r--')
pylab.plot(nthreads_kupe, gnu71_kupe_times, 'g-')
pylab.plot(nthreads_kupe, gnu61_kupe_times, 'g--')
pylab.plot(nthreads_kupe, cray_kupe_times, 'b-')
pylab.plot([1, 40], [2*60+58.957, 2*60+58.957], 'c--')
pylab.legend(['intel kupe', 'intel mahuika', 'gnu 7.1 kupe', 'gnu 6.1 kupe', 'cray kupe', 'cray OpenACC'])
pylab.xlabel('number of OpenMP threads')
pylab.ylabel('execution time [s]')
pylab.title('mahuika upwindCxx -numCells 800 -numSteps 10')
pylab.plot(nthreads_kupe, intel_kupe_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_mahuika, intel_mahuika_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_kupe, gnu71_kupe_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_kupe, gnu61_kupe_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_kupe, cray_kupe_times, 'ko', markerfacecolor='None')

pylab.show()

