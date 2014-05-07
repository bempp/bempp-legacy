ExternalProject_Add(
    Boost
    PREFIX ${EXTERNAL_ROOT}
    # Downloads boost from url -- much faster than svn
    URL http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.bz2/download
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ./bootstrap.sh
    BUILD_COMMAND  ./b2 link=static variant=release --with-test
    INSTALL_COMMAND ./b2 link=static variant=release --with-test --prefix=${EXTERNAL_ROOT} install
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)
# Rerun cmake to capture new boost install
add_recursive_cmake_step(Boost DEPENDEES install)
set(BOOST_ROOT "${EXTERNAL_ROOT}" CACHE INTERNAL "Prefix for Boost install")
# Makes sure those are not in the CACHE, otherwise, new version will not be found
unset(Boost_INCLUDE_DIR CACHE)
unset(Boost_LIBRARY_DIR CACHE)
