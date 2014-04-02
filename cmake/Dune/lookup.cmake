# Needs boost: include this after boost.lookup.
if(Boost_BUILT_AS_EXTERNAL_PROJECT)
  set(depends_on Boost)
endif()
set(default_common_URL 
    http://www.dune-project.org/download/2.2.1/dune-common-2.2.1.tar.gz)
set(default_common_MD5 4001d4c95f06e22ded41abeb063c561c)
set(default_geometry_URL
    http://www.dune-project.org/download/2.2.1/dune-geometry-2.2.1.tar.gz)
set(default_geometry_MD5 35bfc7656549abb2a414ecbd06108384)
set(default_grid_URL
    http://www.dune-project.org/download/2.2.1/dune-grid-2.2.1.tar.gz)
set(default_grid_MD5 21f1a53949c1d21682101552e3b5bc5c)
set(default_localfunctions_URL
    http://www.dune-project.org/download/2.2.1/dune-localfunctions-2.2.1.tar.gz)
set(default_localfunctions_MD5 49a8f85802ff5d9ed917c71181dc2fbd)
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
        URL_HASH MD5=${${prefix}${component}_MD5}
    )
endmacro()

set(cmake_args 
    -DCMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT}
    -DCMAKE_DISABLE_FIND_PACKAGE_MPI=TRUE
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_PROGRAM_PATH:PATH=${EXTERNAL_ROOT}/bin
    -DCMAKE_LIBRARY_PATH:PATH=${EXTERNAL_ROOT}/lib
    -DCMAKE_INCLUDE_PATH:PATH=${EXTERNAL_ROOT}/include
    -DDUNE_USE_ONLY_STATIC_LIBS=TRUE
    -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}
    -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
    -DCMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE}
    -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
    -DCMAKE_C_FLAGS_RELWITHDEBINFO=${CMAKE_C_FLAGS_RELWITHDEBINFO}
    -DCMAKE_CXX_FLAGS_RELWITHDEBINFO=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
)       
set(log_args
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

_get_arguments(common)
ExternalProject_Add(
    dune-common
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS ${depends_on}
    ${common_ARGUMENTS}
    CMAKE_ARGS ${cmake_args}
    ${log_args}
)

_get_arguments(geometry)
ExternalProject_Add(
    dune-geometry
    DEPENDS dune-common
    PREFIX ${EXTERNAL_ROOT}
    ${geometry_ARGUMENTS}
    CMAKE_ARGS ${cmake_args}
        -Ddune-common_DIR=${EXTERNAL_ROOT}/src/dune-common-build
    ${log_args}
)

_get_arguments(grid)
ExternalProject_Add(
    dune-grid
    DEPENDS dune-geometry dune-common
    PREFIX ${EXTERNAL_ROOT}
    ${grid_ARGUMENTS}
    CMAKE_ARGS ${cmake_args}
        -Ddune-common_DIR=${EXTERNAL_ROOT}/src/dune-common-build
        -Ddune-geometry_DIR=${EXTERNAL_ROOT}/src/dune-geometry-build
    ${log_args}
)

_get_arguments(localfunctions)
ExternalProject_Add(
    dune-localfunctions
    DEPENDS dune-grid dune-geometry dune-common
    PREFIX ${EXTERNAL_ROOT}
    PATCH_COMMAND
        ${CMAKE_COMMAND} -DROOT=${EXTERNAL_ROOT}/src/dune-localfunctions
            -P ${CURRENT_LOOKUP_DIRECTORY}/patch-localfunctions.cmake
    ${localfunctions_ARGUMENTS}
    CMAKE_ARGS ${cmake_args}
        -Ddune-common_DIR=${EXTERNAL_ROOT}/src/dune-common-build
        -Ddune-geometry_DIR=${EXTERNAL_ROOT}/src/dune-geometry-build
        -Ddune-grid_DIR=${EXTERNAL_ROOT}/src/dune-grid-build
    ${log_args}
)

ExternalProject_Add(
    dune-foamgrid
    DEPENDS dune-grid dune-geometry dune-common dune-localfunctions
    PREFIX ${EXTERNAL_ROOT}
    GIT_REPOSITORY https://users.dune-project.org/repositories/projects/dune-foamgrid.git 
    PATCH_COMMAND
       ${CMAKE_COMMAND} -E copy_if_different
                        ${CURRENT_LOOKUP_DIRECTORY}/foamgrid-install.cmake
                        ${EXTERNAL_ROOT}/src/dune-foamgrid/CMakeLists.txt
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT}
               -DROOT=${PROJECT_SOURCE_DIR}
    ${log_args}
)
#   URL ${PROJECT_SOURCE_DIR}/contrib/dune/dune-foamgrid

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
  DOWNLOAD_COMMAND dune-common/bin/duneproject < bempp.dune.input
  DEPENDS dune-common dune-foamgrid dune-foamgrid dune-grid dune-geometry dune-localfunctions
  PREFIX ${EXTERNAL_ROOT}
  CMAKE_ARGS -DCMAKE_PROGRAM_PATH:PATH=${EXTERNAL_ROOT}/bin
             -DCMAKE_LIBRARY_PATH:PATH=${EXTERNAL_ROOT}/lib
             -DCMAKE_INCLUDE_PATH:PATH=${EXTERNAL_ROOT}/include
             -DDUNE_USE_ONLY_STATIC_LIBS=TRUE
             -DCMAKE_DISABLE_FIND_PACKAGE_MPI=TRUE
  INSTALL_COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
           ${EXTERNAL_ROOT}/src/dune-bempp-build/config.h
           ${EXTERNAL_ROOT}/include/dune_config.h
)
# Rerun cmake to capture new dune install
add_recursive_cmake_step(dune-bempp
    PACKAGE_NAME Dune
    FOUND_VAR Dune_FOUND
    DEPENDEES install
)
