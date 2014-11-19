# Add a warning if pcre library is not found
find_program(PCRE_CONFIG_EXECUTABLE pcre-config)
if(NOT PCRE_CONFIG_EXECUTABLE)
    if(NOT EXISTS ${EXTERNAL_ROOT}/src/SWIG/pcre-8.36.tar.gz)
        file(WRITE "${EXTERNAL_ROOT}/src/pcre.cmake"
            "file(DOWNLOAD\n"
            "\"ftp://ftp.csx.cam.ac.uk/pub/software/"
                "programming/pcre/pcre-8.36.tar.gz\"\n"
            "\"${EXTERNAL_ROOT}/src/SWIG/pcre-8.36.tar.gz\"\n"
            ")\n"
        )
    endif()
endif()
# Downloads SWIG
ExternalProject_Add(
    SWIG
    PREFIX ${EXTERNAL_ROOT}
    URL http://prdownloads.sourceforge.net/swig/swig-3.0.2.tar.gz
    URL_MD5 62f9b0d010cef36a13a010dc530d0d41
    CONFIGURE_COMMAND ./configure --prefix=${EXTERNAL_ROOT}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    TIMEOUT 10
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)
if(NOT PCRE_CONFIG_EXECUTABLE)
    ExternalProject_Add_Step(SWIG
        pcre
        COMMAND ${CMAKE_COMMAND} -P "${EXTERNAL_ROOT}/src/pcre.cmake"
        COMMAND Tools/pcre-build.sh
        WORKING_DIRECTORY "${EXTERNAL_ROOT}/src/SWIG"
        DEPENDERS configure
        DEPENDEES download
    )

endif()

add_recursive_cmake_step(SWIG DEPENDEES install)

# Makes sure those are not in the CACHE, otherwise, new version will not be found
unset(SWIG_DIR CACHE)
unset(SWIG_EXECUTABLE CACHE)
unset(SWIG_VERSION CACHE)
