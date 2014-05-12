include(PassonVariables)
passon_variables(Armadillo
    FILENAME "${EXTERNAL_ROOT}/src/ArmadilloVariables.cmake"
    PATTERNS
        "CMAKE_[^_]*_R?PATH" "CMAKE_C_.*" "CMAKE_CXX_.*"
        "BLAS_.*" "LAPACK_.*"
    ALSOADD
        "\nset(CMAKE_INSTALL_PREFIX \"${EXTERNAL_ROOT}\" CACHE STRING \"\")\n"
)
if(BLAS_mkl_core_LIBRARY)
    get_filename_component(libdir "${BLAS_mkl_core_LIBRARY}" PATH)
    file(APPEND "${EXTERNAL_ROOT}/src/ArmadilloVariables.cmake"
       "set(CMAKE_LIBRARY_PATH \${CMAKE_LIBRARY_PATH} \"${libdir}\" CACHE PATH \"\")\n"
    )
endif()
ExternalProject_Add(
    Armadillo
    PREFIX ${EXTERNAL_ROOT}
    URL http://sourceforge.net/projects/arma/files/armadillo-4.000.4.tar.gz
    CMAKE_ARGS
        -C "${EXTERNAL_ROOT}/src/ArmadilloVariables.cmake"
        -DCMAKE_BUILD_TYPE=Release
    TIMEOUT 10
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)
# Rerun cmake to capture new armadillo install
add_recursive_cmake_step(Armadillo DEPENDEES install)
