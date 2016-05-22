configure_file(upwindCxx_fitzroy.ll.in upwindCxx_fitzroy.ll)
if (HAVE_OPENMP)
	add_definitions(-DHAVE_OPENMP)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
add_executable(upwindCxx upwind.cxx)
add_executable(upwindAccCxx upwindAcc.cxx)
add_executable(upwindAcc2Cxx upwindAcc2.cxx)

add_test(NAME upwindCxx1
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindCxx1
                     PROPERTIES 
                     ENVIRONMENT "OMP_NUM_THREADS=1"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindCxx2
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindCxx2
                     PROPERTIES 
                     ENVIRONMENT "OMP_NUM_THREADS=2"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindCxx4
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindCxx4
                     PROPERTIES 
                     ENVIRONMENT "OMP_NUM_THREADS=4"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindCxx8
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindCxx8
                     PROPERTIES 
                     ENVIRONMENT "OMP_NUM_THREADS=8"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindCxx
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindCxx
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindAccCxx 
         COMMAND upwindAccCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")

set_tests_properties(upwindAccCxx
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindAcc2Cxx 
         COMMAND upwindAcc2Cxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")

set_tests_properties(upwindAcc2Cxx
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

