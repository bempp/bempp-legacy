# Looks for trilinos. If not found, download and install it.
if(Trilinos_ARGUMENTS)
    cmake_parse_arguments(Trilinos "" "URL;MD5" "" ${Trilinos_ARGUMENTS})
endif()
if(NOT Trilinos_URL)
    set(arguments 
        URL;
        http://trilinos.sandia.gov/download/files/trilinos-11.6.1-Source.tar.bz2
        URL_HASH;
        MD5=b97d882535fd1856599b1c7338f5b45a
    )
elseif(Trilinos_MD5)
    set(arguments URL;${Trilinos_URL};URL_HASH;MD5=${Trilinos_MD5})
else()
    message(FATAL_ERROR "URL specified, but no MD5. Aborting")
endif()

# Create list of dependencies
set(depends_on)
foreach(component TBB Boost SWIG)
  if(${component}_BUILT_AS_EXTERNAL_PROJECT)
    list(APPEND depends_on ${component})
  endif()
endforeach()

set(
  tbb_libraries
  ${TBB_LIBRARY};${TBB_LIBRARY_DEBUG};${TBB_MALLOC_LIBRARY};${TBB_MALLOC_LIBRARY_DEBUG} 
)
# Downloads armadillo
ExternalProject_Add(
    Trilinos
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS ${depends_on}
    ${arguments}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT}
               -DCMAKE_PROGRAM_PATH:PATH=${EXTERNAL_ROOT}/bin
               -DCMAKE_LIBRARY_PATH:PATH=${EXTERNAL_ROOT}/lib
               -DCMAKE_INCLUDE_PATH:PATH=${EXTERNAL_ROOT}/include
               -DBUILD_SHARED_LIBS:BOOL=ON 
               -DTPL_TBB_LIBRARIES:STRING=${tbb_libraries}
               -DTPL_TBB_INCLUDE_DIRS:PATH=${TBB_INCLUDE_DIR}
               -DTPL_ENABLE_BLAS:BOOL=ON
               -DTPL_ENABLE_LAPACK:BOOL=ON
               -DTPL_ENABLE_TBB:BOOL=ON
               -DTPL_ENABLE_Boost:BOOL=ON
               -DTPL_ENABLE_MPI:BOOL=OFF
               -DTPL_Boost_INCLUDE_DIRS:STRING=${Boost_INCLUDE_DIR}
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
               -DTrilinos_ENABLE_PyTrilinos:BOOL=ON
               -DTrilinos_ENABLE_AztecOO:BOOL=ON
               -DTrilinos_ENABLE_Fortran:BOOL=OFF
               -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
               -DRTOp_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
               -DStratimikos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
               -DThyra_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
               -DTpetra_INST_COMPLEX_FLOAT:BOOL=ON
               -DTpetra_INST_FLOAT:BOOL=ON
               -DTPL_ENABLE_MPI:BOOL=${WITH_MPI}
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)
# Rerun cmake to capture new armadillo install
add_recursive_cmake_step(Trilinos Trilinos_FOUND DEPENDEES install)
