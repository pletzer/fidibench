
if (HAVE_OPENMP)
	add_definitions(-DHAVE_OPENMP)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
add_executable(upwind upwind.cxx)
add_test(NAME upwind COMMAND upwind 128 10)
