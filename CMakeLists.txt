cmake_minimum_required (VERSION 2.8)
project(fidibench)

enable_language(Fortran)
enable_testing()

# make sure the default build is of Release type
# use cmake -D CMAKE_BUILD_TYPE=Debug to create a debug version
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(OpenMP)
find_package(PythonExecutable)
find_package(JuliaExecutable)

# Use cmake -D NUM_CELLS=# and -D NUM_TIME_STEPS=# to change default values
set(NUM_CELLS "32" CACHE STRING "Number of cells along each direction when running tests")
set(NUM_TIME_STEPS "10" CACHE STRING "Number of time steps when running tests")
set(NUM_PROCS "8" CACHE STRING "Number of processes when running tests")

if (${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
  # On Fitzroy AIX find_package(MPI) will fail 
  # unless one sets MPI_CXX_INCLUDE_PATH and 
  # MPI_CXX_LIBRARIES
  set(POE "/usr/lpp/ppe.poe")
  if(EXISTS "${POE}/include" AND EXISTS "${POE}/lib")
    set(MPI_INCLUDE_PATH "${POE}/include")
      set(MPI_LIBRARIES "${POE}/lib/libmpi_r.a;${POE}/lib/libvtd_r.a")
      set(MPI_C_INCLUDE_PATH "${MPI_INCLUDE_PATH}")
      set(MPI_C_LIBRARIES "${MPI_LIBRARIES}")
      set(MPI_CXX_INCLUDE_PATH "${MPI_INCLUDE_PATH}")
      set(MPI_CXX_LIBRARIES "${MPI_LIBRARIES}")
      # have not tested that part...
      set(MPI_Fortran_INCLUDE_PATH "${MPI_INCLUDE_PATH}")
      set(MPI_Fortran_LIBRARIES "${MPI_LIBRARIES}")
      message(STATUS "Setting MPI_CXX_INCLUDE_PATH=${MPI_CXX_INCLUDE_PATH}")
      message(STATUS "Setting MPI_CXX_LIBRARIES=${MPI_CXX_LIBRARIES}")
  endif()
endif()

message(STATUS "Looking for MPI")
find_package(MPI)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

# enable profiling and tracing with TAU
if (NOT "${TAU_MAKEFILE}" STREQUAL "")
  if ("${CMAKE_VERSION}" VERSION_LESS 2.8.12)
    get_filename_component(TAU_LIBRARY_DIR "${TAU_MAKEFILE}" PATH)
    get_filename_component(TAU_DIR "${TAU_LIBRARY_DIR}" PATH)
  else ()
    get_filename_component(TAU_LIBRARY_DIR "${TAU_MAKEFILE}" DIRECTORY)
    get_filename_component(TAU_DIR "${TAU_LIBRARY_DIR}" DIRECTORY)
  endif ()
  message("*******************************************************************")
  message("Building with TAU enabled, set the following prior to typing 'make'")
  message("export PATH=${TAU_DIR}/bin:$PATH")
  message("export TAU_MAKEFILE=${TAU_MAKEFILE}")
  set(ENV{PATH} "${TAU_DIR}/bin:$PATH")
  set(ENV{TAU_MAKEFILE} "${TAU_MAKEFILE}")
  message("********************************************************************")
  set(CMAKE_CXX_COMPILER "tau_cxx.sh")
endif()

message(STATUS "==============================================================")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "==============================================================")
message(STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "==============================================================")

add_subdirectory(laplacian)
add_subdirectory(upwind)
