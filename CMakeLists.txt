cmake_minimum_required(VERSION 2.8.11)
project(Bempp CXX)

# Get version string
file(STRINGS "VERSION" Bempp_VERSION LIMIT_COUNT 1)
message("BEM++ Version ${Bempp_VERSION}")

string(REPLACE "." ";" VERSION_LIST ${Bempp_VERSION})
list(GET VERSION_LIST 0 BEMPP_VERSION_MAJOR)
list(GET VERSION_LIST 1 BEMPP_VERSION_MINOR)
list(GET VERSION_LIST 2 BEMPP_VERSION_PATCH)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Needed to find BLAS
enable_language(C)

set(CMAKE_C_FLAGS_RELWITHDEBINFO
    "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
# set(CMAKE_CXX_FLAGS "-Wall -Wnon-virtual-dtor -Wno-sign-compare")

set(CMAKE_C_FLAGS
    "${CMAKE_C_FLAGS} -Wno-deprecated")
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wno-deprecated")

# Tells cmake where to find bem++ specific files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Download and installs the cookoff in the build tree
# From here on we can use cookoff recipes
find_package(GreatCMakeCookOff NO_MODULE PATHS "${PROJECT_SOURCE_DIR}/external/GreatCMakeCookOff" QUIET)
initialize_cookoff()
include(utils)
include(PythonInstall)

# These are the original flags given by the user when compiling Bempp
set(ORIGINAL_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# Makes this a c++14 project
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")


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


add_subdirectory(lib)
add_subdirectory(cython)
add_subdirectory(meshes)
add_subdirectory(doc)
# Installs FindX.cmake files for packages depending on cmake
add_subdirectory(cmake)

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

# Add python files to the installation
install_python(DIRECTORY ${CMAKE_SOURCE_DIR}/python/bempp/api DESTINATION bempp
    FILES_MATCHING PATTERN "*.py")

install_python(FILES ${CMAKE_SOURCE_DIR}/cmake/empty_file DESTINATION bempp
    RENAME __init__.py)

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

add_custom_target(bempp DEPENDS libbempp cython-headers cython)

include(CPackConfig)

