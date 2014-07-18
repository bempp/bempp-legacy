# Needs boost: include this after boost.lookup.
if(TARGET Boost)
  set(depends_on Boost)
endif()
set(default_common_URL
    http://www.dune-project.org/download/2.3.1/dune-common-2.3.1.tar.gz)
set(default_common_SHA256 040cd3811d195631cfa99fab43443d69c1b83f82737b0cd98e4f330ec84051f5)
set(default_geometry_URL 
    http://www.dune-project.org/download/2.3.1/dune-geometry-2.3.1.tar.gz)
set(default_geometry_SHA256 caf8f55b79e3217c3e845a9ada48b51a57f090cbbd4e6994e72067f3449b565c)
set(default_grid_URL
    http://www.dune-project.org/download/2.3.1/dune-grid-2.3.1.tar.gz)
set(default_grid_SHA256 f565d3c2562275cba317adef74f75b0a4f6f130abf4e9e1c34712bc9ab63ab03)
set(default_localfunctions_URL
    http://www.dune-project.org/download/2.3.1/dune-localfunctions-2.3.1.tar.gz)
set(default_localfunctions_SHA256 92c2380f58c7c5f6ff6eb0f4ac694626c3bc81686cbef4534bfb44e351f0b771)

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
        "BLAS_.*" "LAPACK_.*"
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

macro(_get_arguments component)
    set(keyvalues  ${component}_URL;${component}_SHA256)
    cmake_parse_arguments(_ "" "${keyvalues}" "" ${Dune_ARGUMENTS})
    if(__${component}_URL AND NOT __${component}_SHA256)
        message(FATAL_ERROR "${component} given a URL but no SHA256 hash")
    elseif(__${component}_URL AND __${component}_SHA256)
        set(prefix "__")
    else()
        set(prefix "default_")
    endif()
    set(${component}_ARGUMENTS
        URL ${${prefix}${component}_URL}
        URL_HASH SHA256=${${prefix}${component}_SHA256}
    )
endmacro()

_get_arguments(common)
ExternalProject_Add(dune-common
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS ${depends_on}
    ${common_ARGUMENTS}
    CMAKE_ARGS -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
)

_get_arguments(geometry)
ExternalProject_Add(dune-geometry
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS dune-common
    ${geometry_ARGUMENTS}
    CMAKE_ARGS -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
)

include(PatchScript)
set(patchdir "${PROJECT_SOURCE_DIR}/cmake/patches/dune")
create_patch_script(Dune dune_patch_script
    CMDLINE "-p0"
    WORKING_DIRECTORY "${EXTERNAL_ROOT}/src"
    "${patchdir}/grid_yaspgrid.patch"
)

_get_arguments(grid)
ExternalProject_Add(dune-grid
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS dune-geometry dune-common ${grid_depends}
    ${grid_ARGUMENTS}
    PATCH_COMMAND ${dune_patch_script}
    CMAKE_ARGS -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
)

_get_arguments(localfunctions)
ExternalProject_Add(dune-localfunctions
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS dune-geometry dune-common
    ${localfunctions_ARGUMENTS}
    PATCH_COMMAND
        ${CMAKE_COMMAND} -DROOT=${EXTERNAL_ROOT}/src/dune-localfunctions
            -P ${CURRENT_LOOKUP_DIRECTORY}/patch-localfunctions.cmake
    CMAKE_ARGS -C "${EXTERNAL_ROOT}/src/DuneVariables.cmake"
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
)

ExternalProject_Add(dune-foamgrid
    DEPENDS dune-grid dune-geometry dune-common
    PREFIX ${EXTERNAL_ROOT}
    URL ${PROJECT_SOURCE_DIR}/contrib/dune/dune-foamgrid
    PATCH_COMMAND
       ${CMAKE_COMMAND} -E copy_if_different
                        ${CURRENT_LOOKUP_DIRECTORY}/foamgrid-install.cmake
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

# Create fake dune project, with the sole goal of generating a config.h file!
# First, generate a new project
ExternalProject_Add(
    dune-bempp
    DEPENDS dune-foamgrid dune-grid dune-geometry dune-common dune-localfunctions
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
    COMMAND dune-common/bin/duneproject < bempp.dune.input
    WORKING_DIRECTORY ${EXTERNAL_ROOT}/src
    COMMENT Creating fake dune-bempp project
    DEPENDERS configure
)

# Creates a single target for Dune and remove subcomponents from ALL
add_custom_target(Dune ALL
    DEPENDS
        dune-foamgrid dune-grid dune-geometry dune-common dune-localfunctions
        dune-bempp
)
set_target_properties(
    dune-foamgrid dune-grid dune-geometry dune-common dune-localfunctions
    dune-bempp
    PROPERTIES EXCLUDE_FROM_ALL TRUE
)


# Rerun cmake to capture new dune install
add_recursive_cmake_step(dune-bempp
    PACKAGE_NAME Dune
    FOUND_VAR Dune_FOUND
    DEPENDEES install
)
add_dependencies(lookup_dependencies dune-bempp)
