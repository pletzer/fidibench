if (JULIA_EXECUTABLE)
  add_test(NAME upwindJulia 
  	       COMMAND "${JULIA_EXECUTABLE}" 
  	       "${CMAKE_SOURCE_DIR}/upwind/julia/upwind.jl" 
  	       "${NUM_CELLS}" "${NUM_TIME_STEPS}")
  add_test(NAME upwindJulia3
  	       COMMAND "${JULIA_EXECUTABLE}" 
  	       "${CMAKE_SOURCE_DIR}/upwind/julia/upwind3.jl" 
  	       "${NUM_CELLS}" "${NUM_TIME_STEPS}")

  set_tests_properties(upwindJulia upwindJulia3 
  	                   PROPERTIES 
  	                   PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1\\.0|0\\.999]")

endif()
