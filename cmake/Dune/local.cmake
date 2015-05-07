### Figures out duneconfig.h by creating fake dune-bempp
### project.
unset(depends)
foreach(component common geometry grid localfunctions)
    if(TARGET dune-${component})
        list(APPEND depends dune-${component})
    endif()
endforeach()

# This file helps to create a fake dune project
# It answers the question posed by the duneproject script, including the subsidiary
# "this directory already exists..."
file(
  WRITE ${EXTERNAL_ROOT}/src/bempp.dune.input
  "dune-bempp
dune-common dune-geometry dune-grid dune-localfunctions 
1
me@me
y
y
y
y
y
"
)
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

# Create fake dune project, with the sole goal of generating a config.h file!
# First, generate a new project
ExternalProject_Add(
    dune-bempp
    DEPENDS ${depends}
    PREFIX ${EXTERNAL_ROOT}
    DOWNLOAD_COMMAND ""
    CMAKE_ARGS -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD OFF LOG_CONFIGURE ON LOG_BUILD ON
    INSTALL_COMMAND
    /bin/bash -c "cp ${EXTERNAL_ROOT}/src/dune-bempp-build/FC.h ${PROJECT_BINARY_DIR}/include/bempp/common/FC.h && cp ${EXTERNAL_ROOT}/src/dune-bempp-build/config.h ${PROJECT_BINARY_DIR}/include/bempp/common/bempp_dune_config.hpp"
    ${build_args}
)
ExternalProject_Add_Step(dune-bempp
    CREATE_PROJECT
    COMMAND /bin/bash -c "DUNE_CONTROL_PATH=${dune-common_PREFIX} ${DuneProject_PROGRAM} < bempp.dune.input"
    WORKING_DIRECTORY ${EXTERNAL_ROOT}/src
    COMMENT Creating fake dune-bempp project
    DEPENDERS configure
)
set_target_properties(
    dune-bempp
    PROPERTIES EXCLUDE_FROM_ALL TRUE
)

if(TARGET dune-alugrid)
    # Seems to be pb when both dune-bempp and dune-alugrid are run a the same time.
    add_dependencies(dune-bempp dune-alugrid)
    add_dependencies(lookup_dependencies dune-bempp)
elseif(NOT TARGET Dune)
    add_custom_target(Dune ALL DEPENDS dune-bempp)
    add_dependencies(lookup_dependencies Dune)
else()
    add_dependencies(Dune dune-bempp)
endif()
