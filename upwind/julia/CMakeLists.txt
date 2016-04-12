if (JULIA_EXECUTABLE)
  add_test(NAME upwindJulia 
  	       COMMAND "${JULIA_EXECUTABLE}" 
  	       "${CMAKE_SOURCE_DIR}/upwind/julia/upwind.jl" 32 10)
  add_test(NAME upwindJulia3
  	       COMMAND "${JULIA_EXECUTABLE}" 
  	       "${CMAKE_SOURCE_DIR}/upwind/julia/upwind3.jl" 32 10)
endif()
