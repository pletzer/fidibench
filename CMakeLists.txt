cmake_minimum_required (VERSION 2.8)
project(fidibench)

enable_language(Fortran)
enable_testing()

set(OPT_FLAGS "" CACHE STRING "compiler optimization flags")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPT_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPT_FLAGS}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(OpenMP)
find_package(PythonExecutable)
find_package(JuliaExecutable)

# Use cmake -D NUM_CELLS=# and -D NUM_TIME_STEPS=# to change default values
set(NUM_CELLS "32" CACHE STRING "Number of cells along each direction when running tests")
set(NUM_TIME_STEPS "10" CACHE STRING "Number of time steps when running tests")
set(NUM_PROCS "8" CACHE STRING "Number of processes when running tests")

# OpenACC
option(OPENACC "Enable OpenACC support" OFF)
option(INTEL_ANALYZER "Add compiler options for Intel Analyzer" OFF)
option(ALLINEA_MAP "Add compiler options for ARM Allinea map profiler" OFF)

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

  # add optimization flags on AIX
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "XL")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -qarch=pwr6 -qtune=pwr6 -q64 -qextname=get_rstart_fname -qinfo=all -qxflag=oldtab -qfree=f90 -qnoescape -qintlog")
  endif()
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "XL")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -q64 -qarch=pwr6 -qtune=pwr6 -qinfo=all")
  endif()
  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

endif()

message(STATUS "Looking for MPI")
find_package(MPI)
if (MPI_FOUND)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

if (INTEL_ANALYZER)
   message(STATUS "You selected INTEL_ANALYZER")
   message(STATUS "You may need to load the following modules:")
   message(STATUS "ml intel/2017a")
   message(STATUS "ml intel/2017a")
   message(STATUS "ml itac/2017.2.028")
   message(STATUS "source itacvars.sh impi5")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -tcollect $ENV{VT_ADD_LIBS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -tcollect $ENV{VT_ADD_LIBS}")
   set(CMAKE_Fortran_FLAGS 
      "${CMAKE_CXX_FLAGS} -g -tcollect $ENV{VT_ADD_LIBS}")
endif()

if (ALLINEA_MAP)
   message(STATUS "You selected ALLINEA_MAP")
   message(STATUS "You may need to module load forge")
   set(CMAKE_C_FLAGS "-g ${CMAKE_C_FLAGS} -Wl,@/nesi/transit/pletzera/fidibench/allinea-profiler.ld")
   set(CMAKE_CXX_FLAGS "-g ${CMAKE_CXX_FLAGS} -Wl,@/nesi/transit/pletzera/fidibench/allinea-profiler.ld")
   set(CMAKE_Fortran_FLAGS "-g ${CMAKE_Fortran_FLAGS} -Wl,@/nesi/transit/pletzera/fidibench/allinea-profiler.ld")
endif()

if ((NOT OPENACC) AND OPENMP_FOUND)
  add_definitions(-DHAVE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  # There is no OpenMP fortran specific flag, fall back onto C flags
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
else ()
  # OpenACC is incompatible with OpenMP, either one or the other or none
  if (OPENACC)
    message(STATUS "OpenACC enabled")
    # These are the PGI openACC flags -- need to find out what the 
    # corresponding options are for the CRAY, gcc, ... compilers 
    # (if they exist)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -acc -Minfo=acc")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -Minfo=acc")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc -Minfo=acc")
  endif () 
endif ()

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
message(STATUS "C       compiler: ${CMAKE_C_COMPILER}")
message(STATUS "           flags: ${CMAKE_C_FLAGS}")
message(STATUS "==============================================================")
message(STATUS "C++     compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "           flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "==============================================================")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "           flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "==============================================================")

add_subdirectory(scripts)
add_subdirectory(cxx)
add_subdirectory(laplacian)
add_subdirectory(upwind)
