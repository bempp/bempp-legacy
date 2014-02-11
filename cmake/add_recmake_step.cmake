# Sets location where external project are included
if(NOT EXTERNAL_ROOT)
  set(EXTERNAL_ROOT ${CMAKE_BINARY_DIR}/external)
endif(NOT EXTERNAL_ROOT)

include(ExternalProject)

# Adds an external step to an external project to rerun cmake
# This means that an externally installed project will be now catch the newly installed project.
# It expects stuff was installed in EXTERNAL_ROOT
macro(add_rerun_cmake_step name) 
  ExternalProject_Add_Step(
    ${name} reCMake
    COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR} 
                     -DCMAKE_PROGRAM_PATH=${EXTERNAL_ROOT}/bin
                     -DCMAKE_LIBRARY_PATH=${EXTERNAL_ROOT}/lib
                     -DCMAKE_INCLUDE_PATH=${EXTERNAL_ROOT}/include
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    ${ARGN}
  )
endmacro() 

# Avoids anoying cmake warning, by actually using the variables.
# The will be if the appropriate find_* is used. But won't be otherwise.
if(CMAKE_PROGRAM_PATH)
endif()
if(CMAKE_LIBRARY_PATH)
endif()
if(CMAKE_INCLUDE_PATH)
endif()
