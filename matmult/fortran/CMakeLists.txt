if (HAVE_OPENMP)
	add_definitions(-DHAVE_OPENMP)
	message(STATUS "OpenMP enabled in Fortran")
	set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
else()
    message(STATUS "OpenMP disabled in Fortran!")
endif()

add_executable(matmult matmult.F90)

add_test(NAME matmult 
         COMMAND matmult 3 4 5)
set_tests_properties(matmult
                     PROPERTIES
                     ENVIRONMENT "OMP_NUM_THREADS=1"
                     FAIL_REGULAR_EXPRESSION "ERROR")

            
