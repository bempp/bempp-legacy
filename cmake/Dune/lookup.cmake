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

# Needed by dune
#enable_language(Fortran)
macro(find_program_or_fail VARIABLE)
    find_program(${VARIABLE} ${ARGN})
    if(NOT ${VARIABLE})
        message(FATAL_ERROR "Program needed for Dune was not found: ${VARIABLE}")
    endif()
endmacro()
find_package(PkgConfig REQUIRED)
find_program_or_fail(libtoolize_EXECUTABLE
    NAMES libtoolize glibtoolize)
find_program_or_fail(autoconf_EXECUTABLE autoconf)
find_program_or_fail(aclocal_EXECUTABLE aclocal)
find_program_or_fail(automake_EXECUTABLE automake)

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

set(build_args
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)
function(add_to_var VARIABLE PATH)
    include(CMakeParseArguments)
    cmake_parse_arguments(envvar "PREPEND" "" "" ${ARGN})
    if("${${VARIABLE}}" STREQUAL "")
        set(separator "")
    else()
        set(separator ":")
    endif()
    if(envvar_PREPEND)
        set(${${VARIABLE}} "${PATH}${separator}${${VARIABLE}}")
    else()
        set(${${VARIABLE}} "${${VARIABLE}}${separator}${PATH}")
    endif()
endfunction()


macro(_file_args OUTVAR withname)
    unset(result_var)
    foreach(filename ${ARGN})
        add_to_var(result_var "${filename}")
    endforeach()
    if(result_var)
        set(${OUTVAR} "--with-${withname}=\"${result_var}\"")
    endif()
endmacro()

_file_args(blas_args blas ${BLAS_LIBRARIES})
_file_args(lapack_args lapack ${LAPACK_LIBRARIES})

find_program(bash_EXECUTABLE bash REQUIRED)

# Create a script that cmake can call
# This should remove some issues that arises when cmake tries to build
# a complicated command line
function(write_configure_file path)
    get_filename_component(filename "${path}" NAME)
    file(WRITE "${PROJECT_BINARY_DIR}/CMakeFiles/external/${filename}"
        "#!${bash_EXECUTABLE}\n"
        "# Calls configure script for Dune packages\n"
        "export CC=${CMAKE_C_COMPILER}\n"
        "export CFLAGS=\"${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE} "
            "-I${BLAS_INCLUDE_DIR} -I${LAPACK_INCLUDE_DIR}\"\n"
        "export CXX=${CMAKE_CXX_COMPILER}\n"
        "export CXXFLAGS=\"${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}\"\n"
        "export FC=${CMAKE_Fortran_COMPILER}\n"
        "export F77=${CMAKE_Fortran_COMPILER}\n"
        "export F90=${CMAKE_Fortran_COMPILER}\n"
        "export FCFLAGS=\"${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}\"\n"
        "export PKG_CONFIG_PATH=\"$ENV{PKG_CONFIG_PATH}\"\n"
        "\n"
        "./configure"
           " --enable-shared=no"
           " --enable-static=yes"
           " --with-pic"
           " --disable-documentation"
           " --enable-fieldvector-size-is-method"
           " --prefix=\"${EXTERNAL_ROOT}\""
           " --with-blas=\"${BLAS_LIBRARIES}\""
           " --with-blas=\"${LAPACK_LIBRARIES}\""
        ${ARGN}
    )
    file(COPY "${PROJECT_BINARY_DIR}/CMakeFiles/external/${filename}"
        DESTINATION "${EXTERNAL_ROOT}/src/"
        FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
    )
endfunction()

set(configure_command "${EXTERNAL_ROOT}/src/dune_configure.sh")
write_configure_file("${configure_command}")


_get_arguments(common)
ExternalProject_Add(
    dune-common
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS ${depends_on}
    ${common_ARGUMENTS}
    CONFIGURE_COMMAND ${configure_command}
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
    ${build_args}
)

_get_arguments(geometry)
ExternalProject_Add(
    dune-geometry
    DEPENDS dune-common
    PREFIX ${EXTERNAL_ROOT}
    ${geometry_ARGUMENTS}
    CONFIGURE_COMMAND ${configure_command}
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
    ${build_args}
)

find_program(PATCH_EXECUTABLE patch REQUIRED)
_get_arguments(grid)
ExternalProject_Add(
    dune-grid
    DEPENDS dune-geometry dune-common ${grid_depends}
    PREFIX ${EXTERNAL_ROOT}
    ${grid_ARGUMENTS}
    CONFIGURE_COMMAND ${grid_configure_command}
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
    ${build_args}
)
ExternalProject_Add_Step(dune-grid
    PATCH
    COMMAND
        ${PATCH_EXECUTABLE} -p0
            < ${PROJECT_SOURCE_DIR}/cmake/patches/dune/grid_yaspgrid.patch
    COMMAND
        ${PATCH_EXECUTABLE} -p0
            < ${PROJECT_SOURCE_DIR}/cmake/patches/dune/grid_dgfparser.patch
    WORKING_DIRECTORY ${EXTERNAL_ROOT}/src
    DEPENDS
        ${PROJECT_SOURCE_DIR}/cmake/patches/dune/grid_dgfparser.patch
        ${PROJECT_SOURCE_DIR}/cmake/patches/dune/grid_yaspgrid.patch
    DEPENDEES download
    DEPENDERS configure
)

_get_arguments(localfunctions)
ExternalProject_Add(
    dune-localfunctions
    DEPENDS dune-geometry dune-common
    PREFIX ${EXTERNAL_ROOT}
    ${localfunctions_ARGUMENTS}
    CONFIGURE_COMMAND ${configure_command}
    PATCH_COMMAND
        ${CMAKE_COMMAND} -DROOT=${EXTERNAL_ROOT}/src/dune-localfunctions
            -P ${CURRENT_LOOKUP_DIRECTORY}/patch-localfunctions.cmake
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
    ${build_args}
)

ExternalProject_Add(
    dune-foamgrid
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
    CONFIGURE_COMMAND ${configure_command}
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
               ${EXTERNAL_ROOT}/src/dune-bempp/config.h
               ${EXTERNAL_ROOT}/include/dune_config.h
    ${build_args}
)
ExternalProject_Add_Step(dune-bempp
    CREATE_PROJECT
    COMMAND dune-common/bin/duneproject < bempp.dune.input
    COMMAND dune-common/bin/dunecontrol --module=dune-bempp autogen
    WORKING_DIRECTORY ${EXTERNAL_ROOT}/src
    DEPENDERS configure
)

# Rerun cmake to capture new dune install
add_recursive_cmake_step(dune-bempp
    PACKAGE_NAME Dune
    FOUND_VAR Dune_FOUND
    DEPENDEES install
)
add_dependencies(lookup_dependencies dune-bempp)
