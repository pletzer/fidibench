# FiDiBench

FiDiBench is a finite difference suite of codes that can be used to benchmark
hardware performance on HPC and other systems. The code examples are small 
enough to be well understood, typically averaging a few hundred lines of code,
but also relevant to scientific computing, which often involves nearest 
neighbour communication.

FiDiBench can also be used to compare the execution speed obtained by 
implementing a given algorithm in different languages (eg C++ vs Python 
vs Julia).

Example 1: comparing execution speed of C++ vs Fortran (upwind 128 cells/100 time steps)
![alt tag](https://raw.githubusercontent.com/pletzer/fidibench/master/pictures/cxx_vs_fortran.png)

Example 2: comparing OpenACC against OpenMP 
![alt tag](https://raw.githubusercontent.com/pletzer/fidibench/master/pictures/openacc_vs_openmp.png)


## Requirements

* CMake 2.8 or later
* Message Passing Interface, MPI (optional)
* OpenMP (optional)
* C/C++ compiler (optional)
* Python 2.7 with numpy (optional)
* Julia 0.45 or later (optional)

## Build instructions

```bash
cd fidibench
mkdir build
cd build
cmake ..
make
ctest
```

### Running tests

```bash
ctest
```

Some tests will not run if there are missing components, eg Julia scripts will not run if
julia is not installed on your platform.


