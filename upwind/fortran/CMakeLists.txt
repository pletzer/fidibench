add_executable(upwindFortran upwind.f90)
add_test(NAME upwindFortran
         COMMAND upwindFortran
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")

set_tests_properties(upwindFortran
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

