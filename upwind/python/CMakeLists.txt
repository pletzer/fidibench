if (PYTHON_FOUND)
  add_test(NAME upwindPython COMMAND ${PYTHON_EXECUTABLE} upwind.py 128 10)
  add_test(NAME upwindPython3 COMMAND ${PYTHON_EXECUTABLE} upwind3.py 128 10)
endif()
