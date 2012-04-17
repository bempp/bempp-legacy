# Boost
set(BOOST_LIBRARY_DIR "" CACHE PATH "Directory containing Boost Libraries")
set(BOOST_INCLUDE_DIR "" CACHE PATH "Directory containing Boost include files")
include_directories(${BOOST_INCLUDE_DIR})

# BLAS
set(BLAS_LIBRARIES "" CACHE STRING "Semicolon-separated list of full paths to BLAS and LAPACK libs")
#set(ARMADILLO_CMAKE "" CACHE STRING "Directory of Armadillo Cmake file")
find_package( Armadillo REQUIRED REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Threading building blocks
set(TBB_INCLUDE_DIR "" CACHE PATH "Full path to Intel TBB include directory")
set(TBB_LIBRARY_DIR "" CACHE PATH "Full path to the directory with Intel TBB libraries")
find_library(TBB_LIBRARY_DEBUG tbb_debug ${TBB_LIBRARY_DIR})
find_library(TBB_LIBRARY_RELEASE tbb ${TBB_LIBRARY_DIR})
include_directories(${TBB_INCLUDE_DIR})

# Ahmed (optional, used only if WITH_AHMED is set)
set(AHMED_INCLUDE_DIR "" CACHE PATH "Full path to the AHMED include directory")
set(AHMED_LIBRARY "" CACHE FILEPATH "Full path to the AHMED library")
set(METIS_LIBRARY "" CACHE FILEPATH "Full path to the METIS library")

# Trilinos (optional, used only if WITH_TRILINOS is set)
set(TRILINOS_ROOT_DIR "" CACHE PATH "Full path to the Trilinos root directory")
set(HDF5_LIBRARY "" CACHE PATH "Full path to the HDF5 library")
if (WITH_TRILINOS)
    find_library(EPETRA_LIBRARY epetra ${TRILINOS_ROOT_DIR}/lib)
    find_library(AMESOS_LIBRARY amesos ${TRILINOS_ROOT_DIR}/lib)
    find_library(THYRA_CORE_LIBRARY thyracore ${TRILINOS_ROOT_DIR}/lib)
    find_library(STRATIMIKOS_LIBRARY stratimikos ${TRILINOS_ROOT_DIR}/lib)
    find_library(STRATIMIKOS_BELOS_LIBRARY stratimikosbelos ${TRILINOS_ROOT_DIR}/lib)
    find_library(TEUCHOS_LIBRARY teuchos ${TRILINOS_ROOT_DIR}/lib)
    find_library(RTOP_LIBRARY rtop ${TRILINOS_ROOT_DIR}/lib)
    set(TRILINOS_INCLUDE_DIR ${TRILINOS_ROOT_DIR}/include)

    # This should be needed only if Trilinos is compiled with MPI (defines HAVE_MPI)
    set(MPI_INCLUDE_DIR "" CACHE PATH "Full path to the MPI include directory")
    include_directories(${MPI_INCLUDE_DIR})
endif ()

# Dune
find_library(LIB_DUNE_COMMON common ${CMAKE_SOURCE_DIR}/contrib/dune/dune-common/dune/common/.libs)
find_library(LIB_DUNE_GRID grid ${CMAKE_SOURCE_DIR}/contrib/dune/dune-grid/dune/grid/.libs)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-common)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-grid)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-localfunctions)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-foamgrid)

set(ALUGRID_ROOT_DIR "" CACHE PATH "Full path to the ALUGrid root directory")
if (WITH_ALUGRID)
    find_library(ALUGRID_LIBRARY alugrid ${ALUGRID_ROOT_DIR}/lib)
    include_directories(${ALUGRID_ROOT_DIR}/include
                        ${ALUGRID_ROOT_DIR}/include/serial)
endif ()

# OpenCL installation location
if (WITH_OPENCL)
    find_package(OPENCL)
    include_directories(${OPENCL_INCLUDE_DIR})
endif()
