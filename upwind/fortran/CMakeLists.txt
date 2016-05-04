
if (HAVE_OPENMP)
	add_definitions(-DHAVE_OPENMP)
	message(STATUS "OpenMP enabled in Fortran")
	set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
else()
    message(STATUS "OpenMP disabled in Fortran!")
endif()

add_executable(upwindFortran upwind.F95)
add_test(NAME upwindFortran
         COMMAND upwindFortran
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")

set_tests_properties(upwindFortran
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

