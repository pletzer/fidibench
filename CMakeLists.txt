cmake_minimum_required (VERSION 2.8)
project(laplacian)

# make sure the default build is of Release type
# use cmake -D CMAKE_BUILD_TYPE=Debug to create a debug version
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

option(BUILD_WITH_TAU "Build with TAU compilers (only on AIX)" OFF) 

# special treatment for AIX because find_package(MPI) fails on fitzroy
if (${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
  message(STATUS "AIX detected")
  set(CMAKE_CXX_COMPILER "mpCC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -q64")
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -qmaxmem=-1")
  elseif (CMAKE_BUILD_TYPE STREQUAL "Power6")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -qmaxmem=-1 -qtune=pwr6 -qarch=pwr6")
  elseif (CMAKE_BUILD_TYPE STREQUAL "OptimizedPower6")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -qmaxmem=-1 -qtune=pwr6 -qarch=pwr6")
  endif()
else ()
  message(STATUS "Looking for MPI")
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

if (BUILD_WITH_TAU)
  set(CMAKE_CXX_COMPILER "tau_cxx.sh")
endif()

add_executable(testMultiArrayIter testMultiArrayIter.cxx)
add_executable(testCubeDecomp testCubeDecomp.cxx)
add_executable(laplacian laplacian.cxx CmdLineArgParser.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  add_executable(upwindTemplated upwindTemplated.cxx)
endif()

message(STATUS "==============================================================")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "==============================================================")
message(STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "==============================================================")
