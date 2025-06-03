configure_file(upwindF08Gpu.sl.in upwindF08Gpu.sl)
configure_file(upwindF08Gpu_azure.sl.in upwindF08Gpu_azure.sl)

if (PSYCLONE_EXECUTABLE)
    message(STATUS "Running ${PSYCLONE_EXECUTABLE}")
    execute_process(COMMAND ${PSYCLONE_EXECUTABLE} -s psy_omp.py -o ${CMAKE_SRC_DIR}/upwind/fortran/upwindF03Psy.F90 upwindF03.F90
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE psyclone_result
                    OUTPUT_VARIABLE psyclone_output
                    ERROR_VARIABLE psyclone_error)
    message(STATUS "psy_omp.py output: ${psyclone_output}")
    message(STATUS "psy_omp.py error: ${psyclone_error}")
    if (psyclone_result EQUAL 0)
        message(STATUS "Psyclone transformation applied successfully.")
        set(PSYCLONE_TRANSFORMED TRUE)
    endif()
else()
    message(STATUS "Psyclone executable not found. Will not apply OpenMP transformation.")
endif()

if (HAVE_OPENMP)
	add_definitions(-DHAVE_OPENMP)
	message(STATUS "OpenMP enabled in Fortran")
	set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
else()
    message(STATUS "OpenMP disabled in Fortran!")
endif()

add_executable(upwindFortran upwind.F90)
add_executable(upwindF03 upwindF03.F90)
add_executable(upwindF08 upwindF08.F90)
add_executable(upwindF08Gpu upwindF08Gpu.F90)


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

add_test(NAME upwindF08Gpu
         COMMAND upwindF08Gpu
                    "${NUM_CELLS}" "${NUM_TIME_STEPS}")

set_tests_properties(upwindF03 upwindF08 upwindF08Gpu
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

            
