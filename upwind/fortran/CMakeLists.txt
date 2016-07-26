configure_file(upwindFortran_fitzroy.ll.in upwindFortran_fitzroy.ll)
if (HAVE_OPENMP)
	add_definitions(-DHAVE_OPENMP)
	message(STATUS "OpenMP enabled in Fortran")
	set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
else()
    message(STATUS "OpenMP disabled in Fortran!")
endif()

add_executable(upwindFortran upwind.F90)
add_executable(upwindF03 upwindF03.F90)
add_executable(upwindF08 upwindF08.f90)

add_test(NAME upwindFortran1
         COMMAND upwindFortran
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindFortran1
                     PROPERTIES
                     ENVIRONMENT "OMP_NUM_THREADS=1"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindFortran2
         COMMAND upwindFortran
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindFortran2
                     PROPERTIES
                     ENVIRONMENT "OMP_NUM_THREADS=2"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindFortran4
         COMMAND upwindFortran
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindFortran4
                     PROPERTIES
                     ENVIRONMENT "OMP_NUM_THREADS=4"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindFortran8
         COMMAND upwindFortran
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindFortran8
                     PROPERTIES
                     ENVIRONMENT "OMP_NUM_THREADS=8"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")


add_test(NAME upwindF03
         COMMAND upwindF03
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")

add_test(NAME upwindF08
         COMMAND upwindF08
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")

set_tests_properties(upwindF03 upwindF08
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

