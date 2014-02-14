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
macro(add_recursive_cmake_step name check_var) 
  set(cmake_arguments -DCMAKE_PROGRAM_PATH:PATH=${EXTERNAL_ROOT}/bin
                      -DCMAKE_LIBRARY_PATH:PATH=${EXTERNAL_ROOT}/lib
                      -DCMAKE_INCLUDE_PATH:PATH=${EXTERNAL_ROOT}/include
                      -DUSE_OWN_${name}:BOOL=TRUE
                      --no-warn-unused-cli)
  if(NOT "${check_var}" STREQUAL "NOCHECK")
    set(cmake_arguments ${cmake_arguments} -D${name}_REQUIRED:INTERNAL=TRUE)
  endif()
  ExternalProject_Add_Step(
    ${name} reCMake
    COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR} ${cmake_arguments}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    ${ARGN}
  )
  if(${${name}_REQUIRED})
    if(NOT ${${check_var}})
      message(FATAL_ERROR "[${name}] Could not be downloaded and installed")
    endif()
  endif()
  # Sets a variable saying we are building this source externally
  set(${name}_BUILT_AS_EXTERNAL_PROJECT TRUE)
endmacro() 

# Avoids anoying cmake warning, by actually using the variables.
# The will be if the appropriate find_* is used. But won't be otherwise.
if(CMAKE_PROGRAM_PATH)
endif()
if(CMAKE_LIBRARY_PATH)
endif()
if(CMAKE_INCLUDE_PATH)
endif()
