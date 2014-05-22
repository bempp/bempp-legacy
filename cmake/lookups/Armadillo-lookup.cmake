if(Armadillo_ARGUMENTS)
    cmake_parse_arguments(Armadillo
        ""
        "URL;MD5;TIMEOUT"
        ""
        ${Armadillo_ARGUMENTS}
    )
endif()
if(NOT Armadillo_URL)
    set(arguments
        URL;
        http://sourceforge.net/projects/arma/files/armadillo-4.000.4.tar.gz
        URL_HASH;
        MD5=c548089e29ee69e9a6e9bce76d270bea
    )
elseif(Armadillo_MD5)
    set(arguments URL;${Armadillo_URL};URL_HASH;MD5=${Armadillo_MD5})
else()
    message(FATAL_ERROR "URL specified, but no MD5. Aborting")
endif()
if(NOT Armadillo_BUILD_TYPE)
    set(Armadillo_BUILD_TYPE Release)
endif()
if(NOT Armadillo_TIMEOUT)
    set(Armadillo_TIMEOUT 15)
endif()

include(PassonVariables)
passon_variables(Armadillo
    FILENAME "${EXTERNAL_ROOT}/src/ArmadilloVariables.cmake"
    PATTERNS
        "CMAKE_[^_]*_R?PATH" "CMAKE_C_.*" "CMAKE_CXX_.*"
        "BLAS_.*" "LAPACK_.*"
    ALSOADD
        "\nset(CMAKE_INSTALL_PREFIX \"${EXTERNAL_ROOT}\" CACHE STRING \"\")\n"
        "\nset(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} CACHE INTERNAL \"\")\n"
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
    ${arguments}
    TIMEOUT ${Armadillo_TIMEOUT}
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
