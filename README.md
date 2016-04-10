# FiDiBench

FiDiBench is a finite difference suite of codes that can be used to benchmark 
performance on HPC and other systems. The code examples are small enough to be well
understood, typically averaging a few hundred lines of code, but also relevant to 
scientific computing, which often involves nearest neighbour communication. 

FiDiBench can also be used to compare the execution speed obtained by implementing an 
algorithm in different languages (eg C++ vs Python vs Julia).

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




