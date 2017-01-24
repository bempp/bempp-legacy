# Looks for MKL linked to numpy

# Function to compute prefix(es) where libraries are located
# Acceptable prefixes must be a path of the form prefix/lib/library or
# prefix/lib64/library
function(_get_prefixes OUTVAR)
    foreach(library ${ARGN})
        get_filename_component(directory "${library}" PATH)
        get_filename_component(subdir "${directory}" NAME)
        if("${subdir}" STREQUAL "lib" OR "${subdir}" STREQUAL "lib64")
            get_filename_component(directory "${directory}" PATH)
            list(APPEND prefixes "${directory}")
        endif()
    endforeach()
    if(NOT "${prefixes}" STREQUAL "")
        list(REMOVE_DUPLICATES prefixes)
        set(${OUTVAR} ${prefixes} PARENT_SCOPE)
    endif()
endfunction()
function(_get_libdirs OUTVAR)
    foreach(library ${ARGN})
        get_filename_component(directory "${library}" PATH)
        list(APPEND libdirs "${directory}")
    endforeach()
    if(NOT "${libdirs}" STREQUAL "")
        list(REMOVE_DUPLICATES libdirs)
        set(${OUTVAR} ${libdirs} PARENT_SCOPE)
    endif()
endfunction()
# Have to set vendor explicitly, otherwise cmake screws up
function(_get_vendor OUTVAR)
    if("${ARGN}" MATCHES "lp64")
        if("${ARGN}" MATCHES "sequential")
            set(${OUTVAR} Intel10_64lp_seq PARENT_SCOPE)
        else()
            set(${OUTVAR} Intel10_64lp PARENT_SCOPE)
        endif()
    elseif("${ARGN}" MATCHES "mkl_ia32")
        set(${OUTVAR} Intel PARENT_SCOPE)
    elseif("${ARGN}" MATCHES "mkl_em64t")
        set(${OUTVAR} Intel PARENT_SCOPE)
    else()
        set(${OUTVAR} Intel10_32 PARENT_SCOPE)
    endif()
endfunction()

if(NumpyMKL_LIBRARIES)
    set(NumpyMKL_FOUND TRUE)
    set(NUMPYMKL_FOUND TRUE)
    _get_prefixes(NumpyMKL_PREFIXES ${NumpyMKL_LIBRARIES})
    _get_libdirs(NumpyMKL_LIBRARY_DIRS ${NumpyMKL_LIBRARIES})
    _get_vendor(NumpyMKL_VENDOR ${NumpyMKL_LIBRARIES})
    return()
endif()

find_package(PythonInterp)
if(PYTHON_EXECUTABLE)
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} numpy_mkl.py
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(result STREQUAL 0)
        set(NumpyMKL_LIBRARIES ${output})
        _get_prefixes(NumpyMKL_PREFIXES ${NumpyMKL_LIBRARIES})
        _get_libdirs(NumpyMKL_LIBRARY_DIRS ${NumpyMKL_LIBRARIES})
        _get_vendor(NumpyMKL_VENDOR ${NumpyMKL_LIBRARIES})
        # Adds pthread if necessary
        if("${NumpyMKL_LIBRARIES}" MATCHES "mkl"
            AND "${NumpyMKL_LIBRARIES}" MATCHES "thread")
            list(APPEND NumpyMKL_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
        endif()
    endif()
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
    NumpyMKL
    REQUIRED_VARS
        NumpyMKL_LIBRARIES
        NumpyMKL_LIBRARY_DIRS
        NumpyMKL_VENDOR
)
if(NUMPYMKL_FOUND AND NOT NumpyMKL_FOUND)
    set(NumpyMKL_FOUND TRUE)
endif()
set(NumpyMKL_LIBRARIES 
    ${NumpyMKL_LIBRARIES}
    CACHE PATH
    "MKL libraries numpy links to"
)
