if (PYTHON_EXECUTABLE)
  set(ENV{PYTHONPATH} "${CMAKE_SOURCE_DIR}/upwind/python:$ENV{PYTHONPATH}")
  message("PYTHONPATH set to $ENV{PYTHONPATH}")

  add_test(NAME upwindPython 
  	       COMMAND "${PYTHON_EXECUTABLE}" 
  	       "${CMAKE_SOURCE_DIR}/upwind/python/upwind.py" 
  	       "${NUM_CELLS}" "${NUM_TIME_STEPS}")

# window creation fails when using one proc with some implementations
  # add_test(NAME upwindPythonMpi1
  #          COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "1"
  #          "${PYTHON_EXECUTABLE}"
  #          "${CMAKE_SOURCE_DIR}/upwind/python/upwindMPI.py"
  #          "${NUM_CELLS}" "${NUM_TIME_STEPS}")

  add_test(NAME upwindPythonMpiN
           COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "${NUM_PROCS}"
           "${PYTHON_EXECUTABLE}"
           "${CMAKE_SOURCE_DIR}/upwind/python/upwindMPI.py"
           "${NUM_CELLS}" "${NUM_TIME_STEPS}")

  # set_tests_properties(upwindPython upwindPythonMpi1 upwindPythonMpiN
  #                      PROPERTIES 
  #                      PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*(1(\\.000[0-9]*)?|0\\.999[0-9]*)[^0-9.]")
  set_tests_properties(upwindPython upwindPythonMpiN
                       PROPERTIES 
                       PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*(1(\\.000[0-9]*)?|0\\.999[0-9]*)[^0-9.]")
endif()
