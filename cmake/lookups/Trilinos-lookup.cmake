# Looks for trilinos. If not found, download and install it.
if(Trilinos_ARGUMENTS)
    cmake_parse_arguments(Trilinos
        "PYPACKED;PYTHON"
        "URL;MD5;BUILD_TYPE;INSTALL_PREFIX"
        ""
        ${Trilinos_ARGUMENTS}
    )
endif()
if(NOT Trilinos_URL)
    set(arguments
        URL;
        http://trilinos.sandia.gov/download/files/trilinos-11.6.2-Source.tar.bz2
        URL_MD5; 15ea6af5559f872919ff19fe5a322eb6
    )
elseif(Trilinos_MD5)
    set(arguments URL;${Trilinos_URL};URL_MD5;${Trilinos_MD5})
else()
    message(FATAL_ERROR "URL specified, but no MD5. Aborting")
endif()
if(NOT Trilinos_BUILD_TYPE)
    set(Trilinos_BUILD_TYPE Release)
endif()
if(Trilinos_PYTHON)
    set(WITH_PYTHON "ON")
else()
    set(WITH_PYTHON "OFF")
endif()

# Create list of dependencies
set(depends_on)
foreach(component TBB Boost SWIG)
  if(TARGET ${component})
    list(APPEND depends_on ${component})
  endif()
endforeach()

set(tbb_libraries ${TBB_LIBRARY} ${TBB_MALLOC_LIBRARY})

# Create a file of variables which is parsed by cmake when configuring Trilinos.
include(PassonVariables)
passon_variables(Trilinos
    FILENAME "${EXTERNAL_ROOT}/src/TrilinosVariables.cmake"
    PATTERNS
        "CMAKE_[^_]*_R?PATH" "CMAKE_C_.*" "CMAKE_CXX_.*"
        "BLAS_.*" "LAPACK_.*" "SWIG_*"
)
get_filename_component(TPL_TBB_INCLUDE_DIRS "${TBB_INCLUDE_DIR}" PATH)
if(Trilinos_PYPACKED AND Trilinos_PYTHON)
    set(prefix_location "${EXTERNAL_ROOT}/python/PyTrilinos")
    file(APPEND "${EXTERNAL_ROOT}/src/TrilinosVariables.cmake"
        "\nset(PyTrilinos_INSTALL_DIR \"${prefix_location}\"\n"
        "   CACHE PATH \"The path where PyTrilinos will be installed\"\n"
        "   FORCE\n"
        ")\n"
        "\nset(PyTrilinos_INSTALL_PREFIX \"\${PyTrilinos_INSTALL_DIR}\"\n"
        "   CACHE PATH \"The path where PyTrilinos will be installed\"\n"
        "   FORCE\n"
        ")\n"
        "\nset(CMAKE_INSTALL_PREFIX \"\${PyTrilinos_INSTALL_DIR}\" CACHE PATH \"\" FORCE)\n"
    )
else()
    file(APPEND "${EXTERNAL_ROOT}/src/TrilinosVariables.cmake"
        "\nset(Trilinos_INSTALL_INCLUDE_DIR \"include/Trilinos\" CACHE PATH \"\")\n"
        "\nset(CMAKE_INSTALL_PREFIX \"${EXTERNAL_ROOT}\" CACHE PATH \"\")\n"
        "\nlist(APPEND CMAKE_INCLUDE_PATH \"${EXTERNAL_ROOT}/include/trilinos\")\n"
    )
endif()
if(Trilinos_PYTHON)
    passon_variables(Trilinos
        FILENAME "${EXTERNAL_ROOT}/src/TrilinosVariables.cmake"
        PATTERNS "PYTHON_.*"
    )
endif()
file(APPEND "${EXTERNAL_ROOT}/src/TrilinosVariables.cmake"
    "\nset(TPL_TBB_LIBRARIES \"${tbb_libraries}\" CACHE STRING \"\")\n"
    "\nset(TPL_TBB_INCLUDE_DIRS \"${TPL_TBB_INCLUDE_DIRS}\" CACHE STRING \"\")\n"
    "\nset(TPL_BLAS_LIBRARIES \"${BLAS_LIBRARIES}\" CACHE STRING \"\")\n"
    "\nset(TPL_LAPACK_LIBRARIES \"${LAPACK_LIBRARIES}\" CACHE STRING \"\")\n"
    "\nset(TPL_BOOST_INCLUDE_DIRS \"${Boost_INCLUDE_DIR}\" CACHE STRING \"\")\n"
    "\nlist(APPEND CMAKE_PREFIX_PATH \"${EXTERNAL_ROOT}\")\n"
    "\nset(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} CACHE PATH \"\")\n"
)

include(PatchScript)
set(patchdir "${PROJECT_SOURCE_DIR}/cmake/patches/Trilinos")
create_patch_script(Trilinos patch_script
    CMDLINE "-p0"
    WORKING_DIRECTORY "${EXTERNAL_ROOT}/src/Trilinos"
    "${patchdir}/Epetra_ConfigDefs.h.patch"
    "${patchdir}/pytrilinos_eigenvalue_typemap.patch"
)

ExternalProject_Add(
    Trilinos
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS ${depends_on}
    ${arguments}
    CMAKE_ARGS -DBUILD_SHARED_LIBS:BOOL=ON
               -DTPL_ENABLE_BLAS:BOOL=ON
               -DTPL_ENABLE_LAPACK:BOOL=ON
               -DTPL_ENABLE_TBB:BOOL=ON
               -DTPL_ENABLE_Boost:BOOL=ON
               -DTPL_ENABLE_MPI:BOOL=OFF
               -DTrilinos_ENABLE_Amesos:BOOL=ON
               -DTrilinos_ENABLE_Anasazi:BOOL=ON
               -DTrilinos_ENABLE_Belos:BOOL=ON
               -DTrilinos_ENABLE_Epetra:BOOL=ON
               -DTrilinos_ENABLE_EpetraExt:BOOL=ON
               -DTrilinos_ENABLE_RTOp:BOOL=ON
               -DTrilinos_ENABLE_Stratimikos:BOOL=ON
               -DTrilinos_ENABLE_Teuchos:BOOL=ON
               -DTrilinos_ENABLE_Thyra:BOOL=ON
               -DTrilinos_ENABLE_ThyraCore:BOOL=ON
               -DTrilinos_ENABLE_ThyraEpetraAdapters:BOOL=ON
               -DTrilinos_ENABLE_ThyraEpetraExtAdapters:BOOL=ON
               -DTrilinos_ENABLE_ThyraTpetraAdapters:BOOL=ON
               -DTrilinos_ENABLE_Tpetra:BOOL=ON
               -DTrilinos_ENABLE_AztecOO:BOOL=ON
               -DTrilinos_ENABLE_Fortran:BOOL=OFF
               -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
               -DRTOp_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
               -DStratimikos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
               -DThyra_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
               -DTpetra_INST_COMPLEX_FLOAT:BOOL=ON
               -DTpetra_INST_FLOAT:BOOL=ON
               -DTPL_ENABLE_MPI:BOOL=${WITH_MPI}
               -DCMAKE_BUILD_TYPE=${Trilinos_BUILD_TYPE}
               -DTrilinos_ENABLE_PyTrilinos:BOOL=${WITH_PYTHON}
               -C ${EXTERNAL_ROOT}/src/TrilinosVariables.cmake
    PATCH_COMMAND ${patch_script}
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

add_recursive_cmake_step(Trilinos DEPENDEES install)
# If installing in bizarre location (ie all under python package dir), then add
# post-lookup script sot that the package can be found.
if(Trilinos_PYPACKED)
    write_lookup_hook(POST_LOOKUP Trilinos
        "list(FIND CMAKE_PREFIX_PATH \"${prefix_location}\" has_prefix)\n"
        "if(has_prefix EQUAL -1)\n"
        "   list(APPEND CMAKE_PREFIX_PATH \"${prefix_location}\")\n"
        "endif()\n"
    )
endif()
