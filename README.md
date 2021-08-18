# FiDiBench

FiDiBench is a finite difference suite of codes that can be used to benchmark
hardware performance on HPC and other systems. The code examples are small 
enough to be well understood, typically averaging a few hundred lines of code,
but also relevant to scientific computing, which often involves nearest 
neighbour communication.

FiDiBench can also be used to compare the execution speed obtained by 
implementing a given algorithm in different languages (eg C++ vs Python 
vs Julia).

Example 1: comparing C++ compilers (upwind 512 cells, 10 time steps). Note the effect vectorization 
when applying -O3 with the GNU compiler. 
![alt tag](https://raw.githubusercontent.com/pletzer/fidibench/master/pictures/mahuika.png)

Example 2: comparing of C++ vs Fortran (upwind 512 cells/10 time steps)
![alt tag](https://raw.githubusercontent.com/pletzer/fidibench/master/pictures/fortran_vs_c++.png)

Example 3: comparing OpenACC against running on a single CPU 
![alt tag](https://raw.githubusercontent.com/pletzer/fidibench/master/pictures/fidibench_openacc.png)


## Requirements

* CMake 3.1 or later
* Message Passing Interface, MPI (optional)
* OpenMP (optional)
* C/C++ compiler supporting the C++11 standard (optional)
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


