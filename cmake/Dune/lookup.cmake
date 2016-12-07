# Needs boost: include this after boost.lookup.
if(TARGET Boost)
  set(depends_on Boost)
endif()
set(default_common_URL
    http://www.dune-project.org/download/2.4.1/dune-common-2.4.1.tar.gz)
set(default_geometry_URL 
    http://www.dune-project.org/download/2.4.1/dune-geometry-2.4.1.tar.gz)
set(default_grid_URL
    http://www.dune-project.org/download/2.4.1/dune-grid-2.4.1.tar.gz)
set(default_localfunctions_URL
    http://www.dune-project.org/download/2.4.1/dune-localfunctions-2.4.1.tar.gz)


# Create list of library paths.
# They will be added to CMAKE_LIBRARY_PATH
set(library_dirs ${CMAKE_LIBRARY_PATH})
foreach(library ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
   get_filename_component(directory "${library}" PATH)
   list(FIND library_dirs "${directory}" index)
   if(NOT index EQUAL -1)
       list(APPEND library_dirs "${directory}")
   endif()
endforeach()
string(REGEX REPLACE ";" " " library_dirs "${library_dirs}")
string(REGEX REPLACE ";" " " rpath_dirs "${CMAKE_INSTALL_RPATH}")

include(PassonVariables)
passon_variables(Dune
    FILENAME "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    PATTERNS
        "CMAKE_[^_]*_R?PATH" "CMAKE_C_.*" "CMAKE_CXX_.*"
        "BLAS_.*" "LAPACK_.*" "CMAKE_Fortran_.*" "ENABLE_Fortran"
    ALSOADD
        "\nset(CMAKE_INSTALL_PREFIX \"${EXTERNAL_ROOT}\" CACHE STRING \"\")\n"
        "set(CMAKE_LIBRARY_PATH ${library_dirs} \"${EXTERNAL_ROOT}/lib\"\n"
        "    CACHE PATH \"\" FORCE)\n"
        "set(CMAKE_INSTALL_RPATH ${rpath_dirs} CACHE INTERNAL \"\")\n"
        "set(CMAKE_PROGRAM_PATH \"${EXTERNAL_ROOT}/bin\" CACHE PATH \"\")\n"
        "set(BUILD_SHARED_LIBS FALSE CACHE BOOL \"\" FORCE)\n"
        "set(CMAKE_BUILD_TYPE Release CACHE INTERNAL \"\" FORCE)\n"
        "set(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS} -fPIC\" CACHE INTERNAL \"\" FORCE)\n"
        "set(CMAKE_DISABLE_FIND_PACKAGE_HDF5 TRUE CACHE BOOL \"\" FORCE)\n"
        "set(CMAKE_DISABLE_FIND_PACKAGE_MPI TRUE CACHE BOOL \"\" FORCE)\n"
        "set(CMAKE_DISABLE_FIND_PACKAGE_Doxygen TRUE CACHE BOOL \"\" FORCE)\n"
        "if(NOT \"$ENV{FC}\" STREQUAL \"\")\n"
        "    set(ENV{FC} \"$ENV{FC}\")\n"
        "endif()\n"
        "if(NOT \"$ENV{FCFLAGS}\" STREQUAL \"\")\n"
        "    set(ENV{FC}FLAGS \"$ENV{FCFLAGS}\")\n"
        "endif()\n"
        "if(NOT \"$ENV{F90}\" STREQUAL \"\")\n"
        "    set(ENV{F90} \"$ENV{F90}\")\n"
        "endif()\n"
        "if(NOT \"$ENV{F90FLAGS}\" STREQUAL \"\")\n"
        "    set(ENV{F90FLAGS} \"$ENV{F90FLAGS}\")\n"
        "endif()\n"
        "if(NOT \"$ENV{F77}\" STREQUAL \"\")\n"
        "    set(ENV{F77} \"$ENV{F77}\")\n"
        "endif()\n"
        "if(NOT \"$ENV{F77FLAGS}\" STREQUAL \"\")\n"
        "    set(ENV{F77FLAGS} \"$ENV{F77FLAGS}\")\n"
        "endif()\n"
)

macro(_get_arguments component)
    set(keyvalues  ${component}_URL;${component}_MD5)
    cmake_parse_arguments(_ "" "${keyvalues}" "" ${Dune_ARGUMENTS})
    if(__${component}_URL AND NOT __${component}_MD5)
        message(FATAL_ERROR "${component} given a URL but no MD5 hash")
    elseif(__${component}_URL AND __${component}_MD5)
        set(prefix "__")
    else()
        set(prefix "default_")
    endif()
    set(${component}_ARGUMENTS
        URL ${${prefix}${component}_URL}
        URL_MD5 ${${prefix}${component}_MD5}
    )
endmacro()

include(PatchScript)
set(patchdir "${PROJECT_SOURCE_DIR}/cmake/patches/dune")
create_patch_script(Dune-Common dune_common_patch_script
    CMDLINE "-p0"
    WORKING_DIRECTORY "${EXTERNAL_ROOT}/src"
    "${patchdir}/dune_common_cmake.patch"
    "${patchdir}/dune_disable_fortran.patch"
    "${patchdir}/dune_macros.patch"
)

create_patch_script(Dune-Geometry dune_geometry_patch_script
    CMDLINE "-p0"
    WORKING_DIRECTORY "${EXTERNAL_ROOT}/src"
    "${patchdir}/dune_geometry_cmake.patch"
)

create_patch_script(Dune-Grid dune_grid_patch_script
    CMDLINE "-p0"
    WORKING_DIRECTORY "${EXTERNAL_ROOT}/src"
    "${patchdir}/dune_grid_cmake.patch"
)

ExternalProject_Add(dune-common
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS ${depends_on}
    URL http://www.dune-project.org/download/2.4.1/dune-common-2.4.1.tar.gz
    # PATCH_COMMAND ${dune_common_patch_script}
    CMAKE_ARGS -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
)

ExternalProject_Add(dune-geometry
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS dune-common
    URL http://www.dune-project.org/download/2.4.1/dune-geometry-2.4.1.tar.gz
    # PATCH_COMMAND ${dune_geometry_patch_script}
    CMAKE_ARGS -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
)


ExternalProject_Add(dune-grid
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS dune-geometry dune-common ${grid_depends}
    URL http://www.dune-project.org/download/2.4.1/dune-grid-2.4.1.tar.gz
    # PATCH_COMMAND ${dune_grid_patch_script}
    CMAKE_ARGS -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
)

ExternalProject_Add(dune-localfunctions
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS dune-geometry dune-common dune-grid
    URL http://www.dune-project.org/download/2.4.1/dune-localfunctions-2.4.1.tar.gz
    CMAKE_ARGS -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
)


# Creates a single target for Dune and remove subcomponents from ALL
add_custom_target(Dune ALL
    DEPENDS
       dune-grid dune-geometry dune-common dune-localfunctions 
)
set_target_properties(
    dune-grid dune-geometry dune-common dune-localfunctions 
    PROPERTIES EXCLUDE_FROM_ALL TRUE
)


# Rerun cmake to capture new dune install
add_recursive_cmake_step(dune-localfunctions
    PACKAGE_NAME Dune
    FOUND_VAR Dune_FOUND
    DEPENDEES install
)
add_dependencies(lookup_dependencies Dune)
