# Check for lapack and blas
# Look for mkl libraries linked to numpy first and use those if available
find_package(NumpyMKL)
if(NumpyMKL_FOUND)
    if(NumpyMKL_PREFIXES)
        list(INSERT CMAKE_PREFIX_PATH 0 ${NumpyMKL_PREFIXES})
        add_to_rpath(${NumpyMKL_LIBRARIES})
    endif()
    set(BLA_VENDOR ${NumpyMKL_VENDOR})
    list(INSERT CMAKE_LIBRARY_PATH 0 ${NumpyMKL_LIBRARY_DIRS})
    set(BLAS_LIBRARIES ${NumpyMKL_LIBRARIES} CACHE STRING "")
    set(LAPACK_LIBRARIES ${NumpyMKL_LIBRARIES} CACHE STRING "")
    set(LAPACK_FOUND TRUE)
endif()
# Then look for other possible blas libraries. The last few lines should ensure
# we are using the numpy MKLs, if found.
set(BLAS_INCLUDE_DIR_OPTIONAL TRUE)
find_package(CBLAS REQUIRED)
# Check if threading exists in library
include(BlasThreads)
if(NOT NumpyMKL_FOUND)
    # Assume mkls have Lapack. Not sure this is true... but can't use numpy's mkl otherwise.
    find_package(LAPACK REQUIRED)
endif()
if(NOT LAPACK_INCLUDE_DIR)
    find_path(LAPACK_INCLUDE_DIR clapack.h)
endif()
