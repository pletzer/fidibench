
add_executable(testMultiArrayIter testMultiArrayIter.cxx)
add_executable(testCubeDecomp testCubeDecomp.cxx)

if (MPI_FOUND)
  add_executable(laplacian laplacian.cxx CmdLineArgParser.cpp)
  add_test(NAME laplacian2D COMMAND ${MPI_EXEC} -n 8 -numDims 2 -numCells 1024)
  add_test(NAME laplacian3D COMMAND ${MPI_EXEC} -n 8 -numDims 3 -numCells 512)
endif()
