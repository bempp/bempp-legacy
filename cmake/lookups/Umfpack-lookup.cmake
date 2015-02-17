string(REPLACE ";" " " BLAS_LIBS "${BLAS_LIBRARIES}")
string(REPLACE ";" " " LAPACK_LIBS "${LAPACK_LIBRARIES}")

file(WRITE "${PROJECT_BINARY_DIR}/external/src/SuiteSparse_config.mk"
        "CC = ${CMAKE_C_COMPILER}\n"
        "CXX = ${CMAKE_CXX_COMPILER}\n"
        "CF = ${CMAKE_C_FLAGS} \$(TARGET_ARCH) -O3 -fexceptions -fPIC -DNTIMER\n"
        "RANLIB = ranlib\n"
        "ARCHIVE = \$(AR) \$(ARFLAGS)\n"
        "CP = cp -f\n"
        "MV = mv -f\n"
        "F77 = \n"
        "F77FLAGS =\n"
        "F77LIB =\n"
        "LIB = -lm\n"
        "INSTALL_LIB = ${PROJECT_BINARY_DIR}/external/lib\n"
        "INSTALL_INCLUDE = ${PROJECT_BINARY_DIR}/external/include\n"
        "XERBLA =\n"
        "BLAS = ${BLAS_LIBS}\n"
        "LAPACK = ${LAPACK_LIBS}\n"
        "GPU_BLAS_PATH =\n"
        "GPU_CONFIG = \n"
        "UMFPACK_CONFIG = \n"
        "CHOLMOD_CONFIG = \$(GPU_CONFIG) -DNPARTITION\n"
        "SPQR_CONFIG = \n"
        "TBB =\n"
        "CLEAN = *.o *.obj *.ln *.bb *.bbg *.da *.tcov *.gcov gmon.out *.bak *.d *.gcda *.gcno\n")

ExternalProject_Add(
    Umfpack
    PREFIX ${EXTERNAL_ROOT}
    URL http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.3.tar.gz
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND cp ${PROJECT_BINARY_DIR}/external/src/SuiteSparse_config.mk ${PROJECT_BINARY_DIR}/external/src/Umfpack/SuiteSparse_config/
    BUILD_COMMAND make library
    INSTALL_COMMAND make install
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)
  

add_recursive_cmake_step(Umfpack DEPENDEES install)
unset(Umfpack_INCLUDE_DIR CACHE)
unset(Umfpack_LIBRARIES CACHE)
     
