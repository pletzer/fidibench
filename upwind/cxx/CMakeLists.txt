configure_file(upwindCxx_fitzroy.ll.in upwindCxx_fitzroy.ll)
if (HAVE_OPENMP)
	add_definitions(-DHAVE_OPENMP)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
add_executable(upwindCxx upwind.cxx)
add_test(NAME upwindCxx 
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")

set_tests_properties(upwindCxx
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

