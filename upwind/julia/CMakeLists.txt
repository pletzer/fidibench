if (JULIA_EXECUTABLE)
  add_test(NAME upwindJulia 
  	       COMMAND "${JULIA_EXECUTABLE}" 
  	       "${CMAKE_SOURCE_DIR}/upwind/julia/upwind.jl" 
  	       "${NUM_CELLS}" "${NUM_TIME_STEPS}")
  add_test(NAME upwindJulia3
  	       COMMAND "${JULIA_EXECUTABLE}" 
  	       "${CMAKE_SOURCE_DIR}/upwind/julia/upwind3.jl" 
  	       "${NUM_CELLS}" "${NUM_TIME_STEPS}")
endif()
