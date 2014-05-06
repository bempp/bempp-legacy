cmake_minimum_required(VERSION 2.8)
project(bempp CXX)

# Needed to find BLAS
enable_language(C)

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
# set(CMAKE_CXX_FLAGS "-Wall -Wnon-virtual-dtor -Wno-sign-compare")

# Tells cmake where to find bem++ specific files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Download and installs the cookoff in the build tree
# From here on we can use cookoff recipes
include(LookUp-GreatCMakeCookOff)

# Makes this a c++11 project
include(AddCPP11Flags)

# Add rpath information so installed libraries can find each other
include(rpath)

# Macros needed
include(BemppOptions)
include(BemppFindDependencies)

# Tests
if(WITH_TESTS)
    include(PythonVirtualEnv)
    enable_testing()
    add_subdirectory(tests)
    # Adds a bash script to run cpp unit tests
    if(UNIX)
        set(EXECUTABLE "cd ${PROJECT_BINARY_DIR}/tests/unit && ./run_tests")
        configure_file("${GreatCMakeCookOff_DIR}/../scripts/localbash.in.sh"
            "${PROJECT_BINARY_DIR}/CMakeFiles/run_cpptests.sh"
            @ONLY
        )
        file(COPY "${PROJECT_BINARY_DIR}/CMakeFiles/run_cpptests.sh"
            DESTINATION ${PROJECT_BINARY_DIR}
            FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        )
    endif()
endif(WITH_TESTS)

add_subdirectory(lib)
add_subdirectory(python)
add_subdirectory(examples)

# Scratch
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/scratch/CMakeLists.txt")
  add_subdirectory(scratch)
endif()

# Uninstall target
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
