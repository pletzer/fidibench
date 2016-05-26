if (PYTHON_EXECUTABLE)
  set(ENV{PYTHONPATH} "${CMAKE_SOURCE_DIR}/upwind/python:$ENV{PYTHONPATH}")
  message("PYTHONPATH set to $ENV{PYTHONPATH}")

  add_test(NAME upwindPython 
  	       COMMAND "${PYTHON_EXECUTABLE}" 
  	       "${CMAKE_SOURCE_DIR}/upwind/python/upwind.py" 
  	       "${NUM_CELLS}" "${NUM_TIME_STEPS}")

  add_test(NAME upwindPython3 
  	       COMMAND "${PYTHON_EXECUTABLE}" 
  	       "${CMAKE_SOURCE_DIR}/upwind/python/upwind3.py"
  	       "${NUM_CELLS}" "${NUM_TIME_STEPS}")

  add_test(NAME upwindPythonMpi1
           COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "1"
           "${PYTHON_EXECUTABLE}"
           "${CMAKE_SOURCE_DIR}/upwind/python/upwindMPI.py"
           "${NUM_CELLS}" "${NUM_TIME_STEPS}")

  add_test(NAME upwindPythonMpiN
           COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "${MPIEXEC_MAX_NUMPROCS}"
           "${PYTHON_EXECUTABLE}"
           "${CMAKE_SOURCE_DIR}/upwind/python/upwindMPI.py"
           "${NUM_CELLS}" "${NUM_TIME_STEPS}")

  set_tests_properties(upwindPython upwindPython3 upwindPythonMpi1 upwindPythonMpiN
  	                   PROPERTIES 
  	                   PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1\\.0|0\\.999]")
endif()
