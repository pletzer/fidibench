
add_executable(testMultiArrayIter testMultiArrayIter.cxx)
add_executable(testCubeDecomp testCubeDecomp.cxx)

if (MPI_FOUND)
  add_executable(laplacian laplacian.cxx CmdLineArgParser.cpp)
  add_test(NAME laplacian2D COMMAND ${MPIEXEC} -n "${NUM_PROCS}" 
  	                        laplacian -numDims 2 -numCells "${NUM_CELLS}")
  add_test(NAME laplacian3D COMMAND ${MPIEXEC} -n "${NUM_PROCS}"
  	                        laplacian -numDims 3 -numCells "${NUM_CELLS}")
endif()
