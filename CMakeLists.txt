cmake_minimum_required(VERSION 2.8)
project(bempp)
enable_language(Fortran) # for BLAS

# Module Path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Macros needed
include(BemppOptions)
include(BemppFindDependencies)

#add_definitions(-DFIBER_CHECK_ARRAY_BOUNDS)
# set(CMAKE_CXX_FLAGS "-Wall -Wnon-virtual-dtor -Wno-sign-compare")


# Main library
add_subdirectory(lib)

# Examples
add_subdirectory(examples) 

# Scratch
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/scratch/CMakeLists.txt")
  add_subdirectory(scratch)
endif ()

# Tests
if (WITH_TESTS)
  add_subdirectory(tests)
endif (WITH_TESTS)

# Documentation
include(documentation)

# Copy all libraries into the same directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Uninstall target
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall 
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

# Always use the correct rpath

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")
