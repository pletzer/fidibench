
add_executable(testMultiArrayIter testMultiArrayIter.cxx)
add_executable(testCubeDecomp testCubeDecomp.cxx)

if (MPI_FOUND)
  add_executable(laplacian laplacian.cxx CmdLineArgParser.cpp)
endif()

if (OPENMP_FOUND)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  add_executable(upwindCxx upwind.cxx)
endif()
