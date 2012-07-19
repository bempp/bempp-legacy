cmake_minimum_required(VERSION 2.8)
project(bempp)
enable_language(Fortran) # for BLAS

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DFIBER_CHECK_ARRAY_BOUNDS")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DFIBER_CHECK_ARRAY_BOUNDS")

# Module Path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Macros needed
include(BemppOptions)
include(BemppFindDependencies)

# set(CMAKE_CXX_FLAGS "-Wall -Wnon-virtual-dtor -Wno-sign-compare")

############################################################################
### find correct calling convention for dotc

file(READ "${CMAKE_SOURCE_DIR}/cmake/test_zdotc.cpp" ZDOTCTEST)
file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_zdotc/src.cxx" "${ZDOTCTEST}\n")

message(STATUS "Test if ZDOTC works with GFORTRAN calling convention...")

try_run(ZDOT_RESULT ZDOT_RESULT_COMPILED
        ${CMAKE_BINARY_DIR}
	${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_zdotc/src.cxx
	CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${BLAS_LIBRARIES}"
	RUN_OUTPUT_VARIABLE ZDOT_OUTPUT)

if(ZDOT_RESULT_COMPILED)
   message(STATUS "Test Output: ${ZDOT_OUTPUT}")
   if ("${ZDOT_RESULT}" EQUAL 0)
      message(STATUS "GFORTRAN calling convention detected")
   else()
      message(STATUS "GFORTRAN calling convention failed. Test for G77 calling convention")
      file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_zdotc/src.cxx" "${ZDOTCTEST}\n")
      try_run(ZDOTG77_RESULT ZDOTG77_RESULT_COMPILED ${CMAKE_BINARY_DIR}
	      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_zdotc/src.cxx
              CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${BLAS_LIBRARIES}"
	      COMPILE_DEFINITIONS "-DG77CONVENTION"
	      RUN_OUTPUT_VARIABLE ZDOT_OUTPUT)
      message(STATUS "Test Output: ${ZDOT_OUTPUT}")
      if ("${ZDOTG77_RESULT}" EQUAL 0)
         message(STATUS "G77 calling convention detected")
	 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG77CONVENTION")
      else()
         message(WARNING "BEMPP will likely not work in complex mode.")
      endif()
   endif()
else()
   unset(ZDOT_RESULT CACHE)
   unset(ZDOT_RESULT_COMPILED CACHE)
   message(FATAL_ERROR "ZDOTC test does not compile. Check libraries")
endif()

unset(ZDOT_RESULT CACHE)
unset(ZDOT_RESULT_COMPILED CACHE)
unset(ZDOTG77_RESULT CACHE)
unset(ZDOTG77_RESULT_COMPILED CACHE)

############################################################################

# Main library
add_subdirectory(lib)

# Python Support
add_subdirectory(python)

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

# Uninstall target
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall 
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
