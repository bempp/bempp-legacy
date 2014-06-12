ExternalProject_Add(
    Patchelf
    PREFIX ${EXTERNAL_ROOT}
    URL http://releases.nixos.org/patchelf/patchelf-0.8/patchelf-0.8.tar.gz
    URL_HASH SHA256=14af06a2da688d577d64ff8dac065bb8903bbffbe01d30c62df7af9bf4ce72fe
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
