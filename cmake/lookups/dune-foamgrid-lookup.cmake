# Create a script that cmake can call
# This should remove some issues that arises when cmake tries to build
# a complicated command line

if(NOT FOAMGRID_BUILD_TYPE)
    set(FOAMGRID_BUILD_TYPE Release)
endif()


if(NOT TARGET Dune)
    include(PassonVariables)
    passon_variables(Dune
        FILENAME "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
        PATTERNS
            "CMAKE_[^_]*_R?PATH" "CMAKE_C_.*" "CMAKE_CXX_.*"
            "BLAS_.*" "LAPACK_.*" "Dune.*" "dune.*"
        ALSOADD
            "\nset(CMAKE_INSTALL_PREFIX \"${EXTERNAL_ROOT}\" CACHE STRING \"\")\n"
            "set(CMAKE_LIBRARY_PATH ${library_dirs} \"${EXTERNAL_ROOT}/lib\"\n"
            "    CACHE PATH \"\" FORCE)\n"
            "set(CMAKE_INSTALL_RPATH ${rpath_dirs} CACHE INTERNAL \"\")\n"
            "set(CMAKE_PROGRAM_PATH \"${EXTERNAL_ROOT}/bin\" CACHE PATH \"\")\n"
            "set(BUILD_SHARED_LIBS TRUE CACHE BOOL \"\" FORCE)\n"
            "set(CMAKE_BUILD_TYPE Release CACHE INTERNAL \"\" FORCE)\n"
            "set(CMAKE_DISABLE_FIND_PACKAGE_HDF5 TRUE CACHE BOOL \"\" FORCE)\n"
            "set(CMAKE_DISABLE_FIND_PACKAGE_MPI TRUE CACHE BOOL \"\" FORCE)\n"
            "set(CMAKE_DISABLE_FIND_PACKAGE_Doxygen TRUE CACHE BOOL \"\" FORCE)\n"
    )
endif()
unset(depends)
set(depends)
if(TARGET Dune)
    list(APPEND depends Dune)
endif()

ExternalProject_Add(
    dune-foamgrid 
    DEPENDS ${depends}
    PREFIX ${EXTERNAL_ROOT}
    GIT_REPOSITORY https://github.com/bempp/dune-foamgrid.git
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${FOAMGRID_BUILD_TYPE}
               -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
)
add_recursive_cmake_step(dune-foamgrid DEPENDEES install)
