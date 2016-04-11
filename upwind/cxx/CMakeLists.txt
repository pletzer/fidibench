
if (HAVE_OPENMP)
	add_definitions(-DHAVE_OPENMP)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
add_executable(upwindCxx upwind.cxx)
add_test(NAME upwindCxx COMMAND upwind 128 10)
