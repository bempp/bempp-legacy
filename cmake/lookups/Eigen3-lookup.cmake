if(NOT IS_DIRECTORY ${EXTERNAL_ROOT}/include)
    file(MAKE_DIRECTORY ${EXTERNAL_ROOT}/include)
endif()

ExternalProject_Add(Eigen3
    PREFIX ${EXTERNAL_ROOT}
    URL http://bitbucket.org/eigen/eigen/get/3.2.4.tar.bz2
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND cp -R ${EXTERNAL_ROOT}/src/Eigen3/Eigen ${EXTERNAL_ROOT}/include/Eigen
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    )

add_recursive_cmake_step(Eigen3 DEPENDEES install)
unset(Eigen3_INCLUDE_DIR CACHE)
