if(DEFINED BEMPP_MKL_SET_THREADS)
    return()
endif()
file(WRITE "${PROJECT_BINARY_DIR}/CMakeFiles/blas_threads/main.cc"
    "extern \"C\"  void FUNCTION_NAME(int nth);\n"
    "int main(int argc, char const* argv[]) {\n"
    "  FUNCTION_NAME(1);\n"
    "  return 0;\n"
    "};\n"
)
try_compile(BEMPP_MKL_SET_THREADS
    "${PROJECT_BINARY_DIR}/CMakeFiles/blas_threads/intel/"
    "${PROJECT_BINARY_DIR}/CMakeFiles/blas_threads/main.cc"
    COMPILE_DEFINITIONS FUNCTION_NAME=MKL_Set_Num_Threads
    LINK_LIBRARIES ${BLAS_LIBRARIES}
)
try_compile(BEMPP_OPENBLAS_SET_THREADS
    "${PROJECT_BINARY_DIR}/CMakeFiles/blas_threads/intel/"
    "${PROJECT_BINARY_DIR}/CMakeFiles/blas_threads/main.cc"
    COMPILE_DEFINITIONS FUNCTION_NAME=openblas_set_num_threads
    LINK_LIBRARIES ${BLAS_LIBRARIES}
)
try_compile(BEMPP_GOTOBLAS_SET_THREADS
    "${PROJECT_BINARY_DIR}/CMakeFiles/blas_threads/intel/"
    "${PROJECT_BINARY_DIR}/CMakeFiles/blas_threads/main.cc"
    COMPILE_DEFINITIONS FUNCTION_NAME=goto_set_num_threads
    LINK_LIBRARIES ${BLAS_LIBRARIES}
)

unset(WITH_MKL CACHE)
unset(WITH_OPENBLAS CACHE)
unset(WITH_GOTOBLAS CACHE)
if(BEMPP_MKL_SET_THREADS)
    set(WITH_MKL TRUE CACHE BOOL "Using MKL libraries." FORCE)
    message(STATUS "Found MKL_Set_Num_Threads")
elseif(BEMPP_OPENBLAS_SET_THREADS)
    set(WITH_OPENBLAS TRUE CACHE BOOL "Using OpenBLAS libraries." FORCE)
    message(STATUS "Found openblas_set_num_threads")
elseif(BEMPP_GOTOBLAS_SET_THREADS)
    set(WITH_GOTOBLAS TRUE)
    set(WITH_GOTOBLAS TRUE CACHE BOOL "Using GotoBLAS libraries." FORCE)
    message(STATUS "Found goto_set_num_threads")
endif()
