include_directories(${CMAKE_SOURCE_DIR}/cxx/)

add_executable(testMultiArrayIter testMultiArrayIter.cxx)
target_link_libraries(testMultiArrayIter fidibench)
add_executable(testCubeDecomp testCubeDecomp.cxx)
target_link_libraries(testCubeDecomp fidibench)

if (MPI_FOUND)
  add_executable(testStencil2d testStencil2d.cxx)
  target_link_libraries(testStencil2d fidibench)
  add_test(NAME testStencil2d 
  	       COMMAND ${MPIEXEC} -n "${NUM_PROCS}" 
  	            ./testStencil2d -numCells "${NUM_CELLS}")

  add_executable(laplacian laplacian.cxx)
  target_link_libraries(laplacian fidibench)
  add_test(NAME laplacian2D COMMAND ${MPIEXEC} -n "${NUM_PROCS}" 
  	                        ./laplacian -numDims 2 -numCells "${NUM_CELLS}")
  add_test(NAME laplacian3D COMMAND ${MPIEXEC} -n "${NUM_PROCS}"
  	                        ./laplacian -numDims 3 -numCells "${NUM_CELLS}")

  # we need a better test here
  set_tests_properties(laplacian2D laplacian3D
  	                   PROPERTIES 
  	                   PASS_REGULAR_EXPRESSION "[Cc]heck sums:[ ]*")

endif()
