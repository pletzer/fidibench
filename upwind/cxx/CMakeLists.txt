# configure 1, 2, 4, .. thread submission scripts
foreach(i 1 2 4 8 16 32 64)
    configure_file(upwindCxx_fitzroy_${i}.ll.in upwindCxx_fitzroy_${i}.ll)
endforeach(i)
foreach(i 1 2 4 8 16 32)
    configure_file(upwindCxx_pan_${i}.sl.in upwindCxx_pan_${i}.sl)
endforeach(i)
configure_file(upwindMpiCxx_pan_20.sl.in upwindMpiCxx_pan_20.sl)

if (HAVE_OPENMP)
	add_definitions(-DHAVE_OPENMP)
endif()

include_directories(${CMAKE_SOURCE_DIR}/laplacian/cxx/)

set(EXTRA_LIBS "")
if (INTEL_ANALYZER)
  set(EXTRA_LIBS "${EXTRA_LIBS} $ENV{VT_ADD_LIBS}") 
  string(STRIP "${EXTRA_LIBS}" EXTRA_LIBS)
  message(STATUS "Extra libraries are ${EXTRA_LIBS}")
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
add_executable(upwindCxx upwind.cxx
               ${CMAKE_SOURCE_DIR}/laplacian/cxx/CmdLineArgParser.cpp)
target_link_libraries(upwindCxx "${EXTRA_LIBS}")

add_executable(upwindAccCxx upwindAcc.cxx)
target_link_libraries(upwindAccCxx "${EXTRA_LIBS}")

add_executable(upwindAcc2Cxx upwindAcc2.cxx)
target_link_libraries(upwindAcc2Cxx "${EXTRA_LIBS}")

add_executable(upwindMpiCxx upwindMpi.cxx 
               ${CMAKE_SOURCE_DIR}/laplacian/cxx/CmdLineArgParser.cpp
               ${CMAKE_SOURCE_DIR}/laplacian/cxx/Filter.h 
               ${CMAKE_SOURCE_DIR}/laplacian/cxx/CubeDecomp.h  
               ${CMAKE_SOURCE_DIR}/laplacian/cxx/MultiArrayIter.h  
               ${CMAKE_SOURCE_DIR}/laplacian/cxx/writeVTK.h)
target_link_libraries(upwindMpiCxx "${EXTRA_LIBS}")

add_test(NAME upwindCxx1
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindCxx1
                     PROPERTIES 
                     ENVIRONMENT "OMP_NUM_THREADS=1"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindCxx2
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindCxx2
                     PROPERTIES 
                     ENVIRONMENT "OMP_NUM_THREADS=2"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindCxx4
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindCxx4
                     PROPERTIES 
                     ENVIRONMENT "OMP_NUM_THREADS=4"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindCxx8
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindCxx8
                     PROPERTIES 
                     ENVIRONMENT "OMP_NUM_THREADS=8"
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindCxx
         COMMAND upwindCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")
set_tests_properties(upwindCxx
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindAccCxx 
         COMMAND upwindAccCxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")

set_tests_properties(upwindAccCxx
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

add_test(NAME upwindAcc2Cxx 
         COMMAND upwindAcc2Cxx
        "${NUM_CELLS}" "${NUM_TIME_STEPS}")

set_tests_properties(upwindAcc2Cxx
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")

