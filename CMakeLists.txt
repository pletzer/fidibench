cmake_minimum_required (VERSION 2.8)
project(laplacian)

find_package(MPI)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

add_executable(testMultiArrayIter testMultiArrayIter.cxx)
add_executable(testCubeDecomp testCubeDecomp.cxx)
add_executable(laplacian laplacian.cxx CmdLineArgParser.cpp)

