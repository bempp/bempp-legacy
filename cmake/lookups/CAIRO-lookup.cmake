# Downloads Cairo 
ExternalProject_Add(
    CAIRO
    PREFIX ${EXTERNAL_ROOT}
    URL http://cairographics.org/releases/cairo-1.12.16.tar.xz
    URL_HASH SHA256=2505959eb3f1de3e1841023b61585bfd35684b9733c7b6a3643f4f4cbde6d846
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
add_recursive_cmake_step(Cairo DEPENDEES install)

# Makes sure those are not in the CACHE, otherwise, new version will not be found
unset(CAIRO_FOUND CACHE)
unset(CAIRO_INCLUDE_DIRS CACHE)
unset(CAIRO_LIBRARIES CACHE)
