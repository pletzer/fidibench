
add_executable(testMultiArrayIter testMultiArrayIter.cxx)
add_executable(testCubeDecomp testCubeDecomp.cxx)

if (MPI_FOUND)
  add_executable(laplacian laplacian.cxx CmdLineArgParser.cpp)
  add_test(NAME laplacian2D COMMAND ${MPIEXEC} -n 8 laplacian -numDims 2 -numCells 1024)
  add_test(NAME laplacian3D COMMAND ${MPIEXEC} -n 8 laplacian -numDims 3 -numCells 512)
endif()
