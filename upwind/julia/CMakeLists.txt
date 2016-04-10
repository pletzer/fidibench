if (JULIA_FOUND)
  add_test(NAME upwindJulia COMMAND ${JULIA_EXECUTABLE} upwind.jl 128 10)
  add_test(NAME upwindJulia COMMAND ${JULIA_EXECUTABLE} upwind3.jl 128 10)
endif()
