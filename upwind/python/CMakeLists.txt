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
endif()
