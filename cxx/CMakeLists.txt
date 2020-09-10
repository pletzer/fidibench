set(LIB_FILES
  CmdLineArgParser.cpp
  CubeDecomp.cpp
  MultiArrayIter.cpp
  writeVTK.cpp
)

if(MPI_CXX_FOUND)
  set(LIB_FILES ${LIB_FILES} Filter.cpp)
endif()

add_library(fidibench ${LIB_FILES})

