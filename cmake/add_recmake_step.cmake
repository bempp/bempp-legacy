# Sets location where external project are included
if(NOT EXTERNAL_ROOT)
  set(EXTERNAL_ROOT ${CMAKE_BINARY_DIR}/external)
endif(NOT EXTERNAL_ROOT)

include(ExternalProject)

# Adds an external step to an external project to rerun cmake
# This means that an externally installed project will be now catch the newly installed project.
# It expects stuff was installed in EXTERNAL_ROOT
#
# In general, the usage pattern is something like:
#   find_package(something) 
#   if(NOT something_FOUND)
#     ExternalProject_Add(...)
#     add_recursive_cmake_step(something)
#   endif()
#
# This pattern will first attempt to find the package on the system. If it is not found, an external
# project to create it is added, with an extra step to rerun cmake and find the newly installed
# package.
#
# A call to find_package(something) can be protected to make sure that the particular package is
# always downloaded. The pattern is as follows
#
#   if(USE_OWN_something) 
#     find_package(something) 
#   endif()
#   if(NOT something_FOUND)
#     ExternalProject_Add(...)
#     add_recursive_cmake_step(something)
#   endif()
macro(add_recursive_cmake_step name) 
  # If statements limits the depth of the recursive call to cmake
  if(NOT "${PLEASE_NO_RECURSION_ON_${name}}" STREQUAL "PLEASE_NO_RECURSION_ON_${name}")
    ExternalProject_Add_Step(
      ${name} reCMake
      COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR} 
                       -DCMAKE_PROGRAM_PATH:PATH=${EXTERNAL_ROOT}/bin
                       -DCMAKE_LIBRARY_PATH:PATH=${EXTERNAL_ROOT}/lib
                       -DCMAKE_INCLUDE_PATH:PATH=${EXTERNAL_ROOT}/include
                       -DPLEASE_NO_RECURSION_ON_${name}:STRING=PLEASE_NO_RECURSION_ON_${name}
                       -DUSE_OWN_${name}:BOOL=TRUE
                       --no-warn-unused-cli
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      ${ARGN}
    )
  endif()
endmacro() 

# Avoids anoying cmake warning, by actually using the variables.
# The will be if the appropriate find_* is used. But won't be otherwise.
if(CMAKE_PROGRAM_PATH)
endif()
if(CMAKE_LIBRARY_PATH)
endif()
if(CMAKE_INCLUDE_PATH)
endif()
