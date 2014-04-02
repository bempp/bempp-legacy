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
