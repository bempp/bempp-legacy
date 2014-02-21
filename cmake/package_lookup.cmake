# Adds the ability to "look-up" a package
#
# This objective is to define a way to either find a package with find_package and/or,
# depending on choices and circumstances, download and install that package.
#
# The user should only have to include this file and add to their cmake files:
#
#     lookup_package(PACKAGE NAME)
#
# The lookup_package function depends on files in directories in the cmake prefixes paths that have
# the name of the package:
# - ${CMAKE_PREFIX_PATH}/${package}/${package}-lookup.cmake
# - ${CMAKE_PREFIX_PATH}/${package}/${package}-lookup.cmake
# - ${CMAKE_PREFIX_PATH}/${package}/lookup.cmake
#
# These files are included when the function lookup_package is called.
# The files should have the general structure:
#
#     find_package(something) 
#     if(NOT something_FOUND)
#       ExternalProject_Add(...)
#       add_recursive_cmake_step(something)
#     endif()
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
#
#
# This file adds two functions: 
# 
# - lookup_package: calls a lookup file (the ones with the structure above)
# - add_recursive_step: Adds a step to an external project to call cmake after that external project
#   is downloaded/build/installed (to ${EXTERNAL_ROOT}).
#
# It also adds a location in the build directory where external projects are
# downloaded/build/installed. The location is set in the variable EXTERNAL_ROOT.

# Sets location where external project are included
if(NOT EXTERNAL_ROOT)
  set(EXTERNAL_ROOT ${CMAKE_BINARY_DIR}/external)
endif(NOT EXTERNAL_ROOT)
# Makes sure external projects are found by cmake
list(APPEND CMAKE_PREFIX_PATH ${EXTERNAL_ROOT})

include(ExternalProject)

# Looks for a lookup package file and includes it.
function(lookup_package package) 

  foreach(path ${CMAKE_MODULE_PATH})
    list(APPEND cmake_paths ${path}/${package})
  endforeach()
  set(filename ${package}-lookup.cmake)
  find_path(
    LOOKUP_RECIPE ${filename}
    PATHS ${CMAKE_MODULE_PATH} ${cmake_paths} 
    NO_DEFAULT_PATH 
  )

  if("${LOOKUP_RECIPE}" STREQUAL "LOOKUP_RECIPE-NOTFOUND")
    set(filename lookup.cmake)
    find_path(
      LOOKUP_RECIPE lookup.cmake 
      PATHS ${cmake_paths} 
      NO_DEFAULT_PATH 
    )
  endif()
  unset(cmake_paths)


  if("${LOOKUP_RECIPE}" STREQUAL "LOOKUP_RECIPE-NOTFOUND")
    unset(LOOKUP_RECIPE CACHE)
    unset(filename)

    # Checks if package is required
    list(FIND ARGN REQUIRED is_required)
    if(${is_required} LESS 0)
      message(FATAL_ERROR "Could not find recipe to lookup ${package} -- ${LOOKUP_RECIPE}")
    endif()
    message(STATUS "Could not find recipe to lookup ${package} -- ${LOOKUP_RECIPE}")


    return()

  endif()

  set(CURRENT_LOOKUP_DIRECTORY ${LOOKUP_RECIPE})
  include(${LOOKUP_RECIPE}/${filename})
  unset(CURRENT_LOOKUP_DIRECTORY)
  unset(filename)
  unset(LOOKUP_RECIPE CACHE)

endfunction()


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
