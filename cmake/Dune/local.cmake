### Adds dune-foamgrid and figures out duneconfig.h by creating fake dune-bempp
### project.
if(NOT Dune_FOUND)
    # Can only be run if Dune is already found.
    return()
endif()
unset(depends)
foreach(component common geometry grid localfunctions)
    if(TARGET dune-${component})
        list(APPEND depends dune-${component})
    endif()
endforeach()
ExternalProject_Add(dune-foamgrid
    DEPENDS ${depends}
    PREFIX ${EXTERNAL_ROOT}
    URL ${PROJECT_SOURCE_DIR}/contrib/dune/dune-foamgrid
    PATCH_COMMAND
       ${CMAKE_COMMAND} -E copy_if_different
                        ${CMAKE_CURRENT_LIST_DIR}/foamgrid-install.cmake
                        ${EXTERNAL_ROOT}/src/dune-foamgrid/CMakeLists.txt
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT}
)

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
    DEPENDS dune-foamgrid
    PREFIX ${EXTERNAL_ROOT}
    DOWNLOAD_COMMAND ""
    CMAKE_ARGS -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD OFF LOG_CONFIGURE ON LOG_BUILD ON
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
               ${EXTERNAL_ROOT}/src/dune-bempp-build/config.h
               ${EXTERNAL_ROOT}/include/dune_config.h
    ${build_args}
)
ExternalProject_Add_Step(dune-bempp
    CREATE_PROJECT
    COMMAND ${DuneProject_PROGRAM} < bempp.dune.input
    WORKING_DIRECTORY ${EXTERNAL_ROOT}/src
    COMMENT Creating fake dune-bempp project
    DEPENDERS configure
)
set_target_properties(
    dune-foamgrid dune-bempp
    PROPERTIES EXCLUDE_FROM_ALL TRUE
)

if(NOT TARGET Dune)
    add_custom_target(Dune ALL DEPENDS dune-bempp dune-foamgrid)
    add_dependencies(lookup_dependencies Dune)
else()
    add_dependencies(Dune dune-foamgrid dune-bempp)
endif()
