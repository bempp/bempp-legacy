
if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    file(WRITE "${EXTERNAL_ROOT}/src/ffc_install.sh" 
        "cd ${EXTERNAL_ROOT}/src/FFC\n" 
        "MACOSX_DEPLOYMENT_TARGET=10.9 ${PROJECT_BINARY_DIR}/localpython.sh ./setup.py install --prefix ${EXTERNAL_ROOT} --install-lib ${EXTERNAL_ROOT}/python\n")
elseif("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
    file(WRITE "${EXTERNAL_ROOT}/src/ffc_install.sh" 
        "cd ${EXTERNAL_ROOT}/src/FFC\n" 
        "${PROJECT_BINARY_DIR}/localpython.sh ./setup.py install --prefix ${EXTERNAL_ROOT} --install-lib ${EXTERNAL_ROOT}/python\n")
endif()

unset(depends)
set(depends)
if(TARGET SWIG)
    list(APPEND depends SWIG)
endif()


ExternalProject_Add(FFC
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS ${depends}
    URL https://bitbucket.org/fenics-project/ffc/downloads/ffc-1.5.0.tar.gz
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND sh "${EXTERNAL_ROOT}/src/ffc_install.sh"
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)

add_recursive_cmake_step(FFC DEPENDEES install)


