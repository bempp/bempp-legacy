ExternalProject_Add(
    Patchelf
    PREFIX ${EXTERNAL_ROOT}
    URL http://releases.nixos.org/patchelf/patchelf-0.8/patchelf-0.8.tar.gz
    URL_MD5 407b229e6a681ffb0e2cdd5915cb2d01
    CONFIGURE_COMMAND ./configure --prefix=${EXTERNAL_ROOT}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    TIMEOUT 20
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)
