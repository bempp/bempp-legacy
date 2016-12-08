# First attempts using pkg-config. It is not unlikely to work and should be
# quite reliable on linuxes.

# Declares the variables BLAS_INCLUDE_DIR and BLAS_INCLUDE_FILENAME
# The latter could be cblas.h, mkl.h, or Accelerate.h. It can be used to
# include the right header: "#include BLAS_INCLUDE_FILENAME"

#####
#####  Macros corresponding to each step are defined first
#####  And the actual code exectution is below the macros
#####
macro(found_or_not_found_blas)
    if(NOT CBLAS_FIND_QUIETLY)
        if(BLAS_LIBRARIES)
            message(STATUS "Found blas libraries ${BLAS_LIBRARIES}")
        endif()
        if(BLAS_INCLUDE_DIR)
            message(STATUS "Found blas include ${BLAS_INCLUDE_DIR}")
        endif()
    endif()
    if(CBLAS_FIND_REQUIRED)
        if(NOT BLAS_LIBRARIES AND NOT BLAS_INCLUDE_DIR)
            message(FATAL_ERROR "Could not find a blas library")
        elseif(NOT BLAS_LIBRARIES)
            message(FATAL_ERROR "Could not figure out blas library")
        elseif(NOT BLAS_INCLUDE_DIR AND NOT BLAS_INCLUDE_DIR_OPTIONAL)
            message(FATAL_ERROR "Could not figure out blas include dir")
        endif()
    endif()
endmacro()

function(_look_for_blas_libraries)
    include(FindPkgConfig)
    pkg_search_module(CBLAS cblas)
    if(NOT CBLAS_FOUND AND PKG_CONFIG_FOUND)
        pkg_search_module(ATLAS atlas)
    else()
      set(BLAS_LIBRARIES ${CBLAS_LIBRARIES})
      set(BLAS_INCLUDE_DIR ${CBLAS_INCLUDE_DIRS})
    endif()

    # include -pthread so MKL can be included on Ubuntu + enthought
    if(use_pthread_flag)
      set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
      set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -pthread")
    endif()
    find_package(BLAS QUIET)
    if(use_pthread_flag)
      set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
      set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -pthread")
    endif()

    # Figures out atlas if necessary
    if(BLAS_atlas_LIBRARY)
        get_filename_component(atlas_dir "${BLAS_atlas_LIBRARY}" PATH)
        get_filename_component(atlas_ext "${BLAS_atlas_LIBRARY}" EXT)
        find_library(BLAS_atlas_cblas_LIBRARY NAMES libcblas${atlas_ext}
            HINTS "${atlas_dir}"
        )
        if(BLAS_atlas_cblas_LIBRARY)
          set(BLAS_FOUND TRUE)
            set(BLAS_LIBRARIES
                "${BLAS_atlas_cblas_LIBRARY}"
                "${BLAS_atlas_LIBRARY}"
            )
        endif()
    endif()

    # Try open-blas
    if(NOT BLAS_LIBRARIES OR NOT BLAS_INCLUDE_DIR)
        find_package(OpenBLAS QUIET
          PATHS /usr/local $ENV{OpenBLAS_HOME} $ENV{OpenBLAS}
          PATHS_SUFFIXES openblas openblas-base
        )
        if(OpenBLAS_FOUND)
            if(NOT EXISTS "${OpenBLAS_LIBRARIES}")
              string(REPLACE "'" "" OpenBLAS_LIBRARIES ${OpenBLAS_LIBRARIES})
            endif()
            if(NOT EXISTS "${OpenBLAS_INCLUDE_DIRS}")
              string(REPLACE "'" "" OpenBLAS_INCLUDE_DIRS ${OpenBLAS_INCLUDE_DIRS})
            endif()
            set(BLAS_LIBRARIES ${OpenBLAS_LIBRARIES})
            set(BLAS_INCLUDE_DIR ${OpenBLAS_INCLUDE_DIRS})
        endif()
    endif()
    if(BLAS_FOUND)
      set(BLAS_FOUND TRUE PARENT_SCOPE)
    endif()
    if(BLAS_LIBRARIES)
      set(BLAS_LIBRARIES "${BLAS_LIBRARIES}" PARENT_SCOPE)
    endif()
    if(BLAS_INCLUDE_DIR)
      set(BLAS_INCLUDE_DIR "${BLAS_INCLUDE_DIR}" PARENT_SCOPE)
    endif()
endfunction()

function(_look_for_include_directories)
    # Adds BLAS_INCLUDE_DIR
    set(directories)
    foreach(path ${ARGN})
        string(REGEX REPLACE 
            "(.*)/lib(64)?/.*" "\\1/include" current "${path}")
        if(NOT "${current}" STREQUAL "/include"
            AND IS_DIRECTORY "${current}")
            list(APPEND directories "${current}")
        endif()
    endforeach()
    if(NOT "${directories}" STREQUAL "")
        list(REMOVE_DUPLICATES directories)
        set(${OUTVAR} ${directories} PARENT_SCOPE)
    endif()
    # find_package blas does not look for cblas.h
    find_path(BLAS_INCLUDE_DIR
        NAMES cblas.h mkl.h Accelerate.h
        HINTS ${directories}
        PATH_SUFFIXES Headers
    )
    set(BLAS_INCLUDE_DIR ${BLAS_INCLUDE_DIR} PARENT_SCOPE)
endfunction()


macro(_slot_pthread_into_mkl)
    if(use_pthread_flag AND "${BLAS_LIBRARIES}" MATCHES "mkl")
        if(NOT CMAKE_LIB_THREADS_INIT)
            find_library(PTHREAD_LIBRARY pthread)
            if(PTHREAD_LIBRARY)
                set(CMAKE_THREAD_LIBS_INIT ${PTHREAD_LIBRARY})
            endif()
        endif()
        # Add -pthread to CMAKE_C(XX)_FLAGS
        set(BLAS_CMAKE_C_FLAGS -pthread)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
        # Add CMAKE_THREAD_LIBS_INIT right before iomp5
        set(index 0)
        foreach(library ${BLAS_LIBRARIES})
            if("${library}" MATCHES "iomp5")
                list(INSERT BLAS_LIBRARIES ${index} "${CMAKE_THREAD_LIBS_INIT}")
                break()
            elseif("${library}" STREQUAL "${CMAKE_THREAD_LIBS_INIT}")
                break()
            endif()
            math(EXPR index "${index} + 1")
        endforeach()
    endif()
endmacro()

#####
#####  The following code is executed when this file is included
#####
if(NOT DEFINED use_pthread_flag)
    include(CheckCCompilerFlag)
    check_c_compiler_flag(-pthread use_pthread_flag)
    set(use_pthread_flag ${use_pthread_flag} CACHE INTERNAL "")
endif()

if(NOT BLAS_LIBRARIES)
  _look_for_blas_libraries()
endif()

if(BLAS_LIBRARIES)
    _slot_pthread_into_mkl()
    if(NOT BLAS_INCLUDE_DIR)
        _look_for_include_directories(${BLAS_LIBRARIES})
    endif()
endif()

found_or_not_found_blas()

if(EXISTS "${BLAS_INCLUDE_DIR}/cblas.h")
    set(BLAS_INCLUDE_FILENAME cblas.h)
elseif(EXISTS "${BLAS_INCLUDE_DIR}/mkl.h")
    set(BLAS_INCLUDE_FILENAME mkl.h)
elseif(APPLE AND EXISTS "${BLAS_INCLUDE_DIR}/Accelerate.h")
    set(BLAS_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -framework accelerate")
    set(BLAS_INCLUDE_FILENAME Accelerate.h)
endif()
