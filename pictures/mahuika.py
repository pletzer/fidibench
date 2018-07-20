import matplotlib
font = {'family' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
from matplotlib import pylab

# for nth in 36 24 18 12 4 2 1; do export OMP_NUM_THREADS=$nth; time ./upwind/cxx/upwindCxx -numCells 512 -numSteps 10 ; done

# upwindCxx -numCells 512 -numSteps 10
nthreads_mahuika = [1, 2, 4, 12, 18, 24, 36]

#
# mahuika
#

# FC=mpiifort CXX=mpiicc CC=mpiicc cmake CMAKE_CXX_FLAGS="-O3 -ipo -xHost" CMAKE_Fortran_FLAGS="-O3 -ipo -xHost" ../..
intel_2017a_O3_xHost_times = [2*60+31.058, 4*60+7.013, 3*60+7.207, 1*60+24.890, 1*60+6.686, 46.987, 35.067]

# intel default optimization -O2
# FC=mpiifort CXX=mpiicc CC=mpiicc cmake CMAKE_CXX_FLAGS="-O2" CMAKE_Fortran_FLAGS="-O2" ../..
intel_2017a_O2_times = [2*60+31.058, 4*60+7.013, 3*60+7.207, 1*60+24.890, 1*60+6.686, 46.987, 35.067]


# FC=mpif90 CXX=mpicxx CC=mpicc cmake CMAKE_CXX_FLAGS="-O3 -ffast-math -funroll-loops -march=native" CMAKE_Fortran_FLAGS="-O3 -ffast-math -funroll-loops -march=native" ../..
gimkl_O3_fastmath_times = [ 10*60+35.306, 6*60+50.439, 3*60+50.842, 1*60+39.565, 1*60+10.050, 56.322, 40.158]

# Cray compilers
# module load PrgEnv-cray slurm
# FC=ftn CXX=CC CC=cc cmake CMAKE_CXX_FLAGS="-O3 -hfp3 -homp" CMAKE_Fortran_FLAGS="-O3 -hfp3 -homp" ../..
# NEED TO RERUN WHEN SLURM IS WORKING
cray_O3_hpf3_times = [2*60+27.971, 2*60+36.865, 2*60+2.777, 2*60+4.060, 1*60+48.825, 2*60+0.174, 1*60+48.160]

pylab.figure(1)
pylab.plot(nthreads_mahuika, intel_2017a_O3_xHost_times, 'r-')
pylab.plot(nthreads_mahuika, gimkl_O3_fastmath_times, 'g-')
pylab.plot(nthreads_mahuika, cray_O3_hpf3_times, 'b-')
pylab.legend(['intel', 'gnu', 'cray'])
pylab.xlabel('number of OpenMP threads')
pylab.ylabel('execution time [s]')
pylab.title('mahuika upwindCxx -numCells 512 -numSteps 10')
pylab.plot(nthreads_mahuika, intel_2017a_O3_xHost_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_mahuika, gimkl_O3_fastmath_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_mahuika, cray_O3_hpf3_times, 'ko', markerfacecolor='None')

pylab.show()


# FC=mpif90 CXX=mpicxx CC=mpicc cmake CMAKE_CXX_FLAGS="-O2" CMAKE_Fortran_FLAGS="-O2" ../..
gimkl_O2_times = [ 10*60+35.237, 6*60+57.485, 3*60+42.968, 1*60+38.140, 1*60+10.943, 53.456, 43.934]


# Cray compilers 
# module load PrgEnv-cray slurm
# FC=ftn CXX=CC CC=cc cmake CMAKE_CXX_FLAGS="-O3 -hfp3 -homp" CMAKE_Fortran_FLAGS="-O3 -hfp3 -homp" ../..
# NEED TO RERUN WHEN SLURM IS WORKING
cray_O3_hpf3_times = [2*60+27.971, 2*60+36.865, 2*60+2.777, 2*60+4.060, 1*60+48.825, 2*60+0.174, 1*60+48.160]
#
# kupe
#
# module load craype-x86-skylake
# FC=ftn CXX=CC CC=cc cmake CMAKE_CXX_FLAGS="-O3 -hfp3" CMAKE_Fortran_FLAGS="-O3 -hfp3" ../..
nthreads_kupe = [1, 2, 4, 10, 20, 40]
cray_times_skylake = [2*60+12.263, 1*60+18.687, 42.272, 20.423, 13.949, 8.860]
cray_times_broadwell = [2*60+14.858, 1*60+18.482, 42.047, 20.360, 13.960, 8.649]

# FC="ftn -O3 -hfp3 -h preferred_vector_width=128" CXX="CC -O3 -hfp3 -h preferred_vector_width=128" CC=cc cmake ../..
cray_vec128_times = [2*60+16.425, 1*60+20.486, 43.117 , 20.741 , 14.675, 8.808]
cray_vec252_times = [2*60+14.482, 1*60+20.850, 43.269, 20.902, 14.520, 8.783]
cray_vec512_times = [2*60+16.082, 1*60+20.627, 43.121, 20.796, 14.589, 8.882]

# gcc  FC=ftn CXX=CC cmake CMAKE_CXX_FLAGS="-O3 -ffast-math -funroll-loops -march=native" CMAKE_Fortran_FLAGS="-O3 -ffast-math -funroll-loops -march=native" ../..
gcc_4_9_3_times = [9*60+27.759, 4*60+58.405, 2*60+39.420, 1*60+8.620, 41.452, 23.518]
gcc_5_3_0_times = [9*60+20.803, 4*60+51.038, 2*60+35.920, 1*60+7.165, 39.825, 22.806]
gcc_6_1_0_times = [9*60+32.201, 7*60+50.232, 2*60+38.858, 3*60+16.282, 1*60+42.172, 1*60+32.811]
gcc_7_1_0_times = [9*60+28.795, 8*60+2.616, 2*60+37.504, 3*60+13.653, 1*60+41.648, 1*60+32.245]

# intel FC=ftn CXX=CC CC=cc cmake -DMPI_CXX_COMPILER=CC -DMPI_Fortran_COMPILER=ftn -DCMAKE_CXX_FLAGS="-O3 -ipo -xHost" -DCMAKE_Fortran_FLAGS="-O3 -ipo -xHost" ../..
intel_kupe_times = [2*60+1.771, 3*60+4.172, 2*60+15.646, 59.403, 54.169, 35.918]

pylab.figure(2)
pylab.plot(nthreads_mahuika, intel_2017a_O3_xHost_times, 'r-', linewidth=2)
pylab.plot(nthreads_kupe, intel_kupe_times, 'r-')
pylab.plot(nthreads_kupe, gcc_4_9_3_times, 'g-')
pylab.plot(nthreads_kupe, gcc_5_3_0_times, 'g--')
pylab.plot(nthreads_kupe, gcc_6_1_0_times, 'g-.')
pylab.plot(nthreads_kupe, gcc_6_1_0_times, 'g:')
pylab.plot(nthreads_kupe, cray_times_skylake, 'b-')
pylab.legend(['mahuika-intel', 'kupe-intel', 'kupe gcc 4.9.3', 'kupe gcc 5.3.0', 'kupe gcc 6.1.0', 'kupe gcc 7.1.0', 'kupe-cray'])
pylab.plot(nthreads_mahuika, intel_2017a_O3_xHost_times, 'r-', linewidth=2)
pylab.plot(nthreads_kupe, intel_kupe_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_kupe, gcc_4_9_3_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_kupe, gcc_5_3_0_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_kupe, gcc_6_1_0_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_kupe, gcc_6_1_0_times, 'ko', markerfacecolor='None')
pylab.plot(nthreads_kupe, cray_times_skylake, 'ko', markerfacecolor='None')
pylab.xlabel('number of OpenMP threads')
pylab.ylabel('execution time [s]')
pylab.title('mahuika vs kupe: upwindCxx -numCells 512 -numSteps 10')

pylab.show()


# hyperthreading
# getting 31.231


