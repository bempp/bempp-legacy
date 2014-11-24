cmake_minimum_required(VERSION 2.8.6)
project(Bempp CXX)
set(Bempp_VERSION 2.0)

# Needed to find BLAS
enable_language(C)

set(CMAKE_C_FLAGS_RELWITHDEBINFO
    "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
# set(CMAKE_CXX_FLAGS "-Wall -Wnon-virtual-dtor -Wno-sign-compare")

# Tells cmake where to find bem++ specific files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Download and installs the cookoff in the build tree
# From here on we can use cookoff recipes
include(LookUp-GreatCMakeCookOff)

# Makes this a c++11 project
include(cpp11_project)


# Macros to install python
include(PythonInstall)
# Defines where to install libraries, headers
include(InstallPaths)
# Add rpath information so installed libraries can find each other
include(rpath)
add_to_rpath("${LIBRARY_INSTALL_PATH}")

# Macros needed
include(BemppOptions)
include(BemppFindDependencies)

# Documentation target
if(DOXYGEN_FOUND OR SPHINX_FOUND)
    add_custom_target(documentation)
endif()


# Tests
if(WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif(WITH_TESTS)

add_subdirectory(lib)
add_subdirectory(python)
add_subdirectory(examples)
add_subdirectory(meshes)

# Uninstall target
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/uninstall.in.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} 
        -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)

# Exports all Bempp so other packages can access it
include(BemppExport)

# An install script to brutally modify the rpaths
# It ensures that python extensions and libraries can see each other
# Importantly, it works for external projects as well as bemmp objects
configure_file("${PROJECT_SOURCE_DIR}/cmake/PostInstall_FixRpaths.in.cmake"
    "${PROJECT_BINARY_DIR}/post_install/FixRPaths.cmake"
    @ONLY
)
file(WRITE "${PROJECT_BINARY_DIR}/post_install/CMakeLists.txt"
    "install(SCRIPT FixRPaths.cmake)"
)
add_subdirectory("${PROJECT_BINARY_DIR}/post_install"
    "${PROJECT_BINARY_DIR}/post_install-build")
