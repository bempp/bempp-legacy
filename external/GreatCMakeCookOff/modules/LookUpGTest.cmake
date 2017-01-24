# Installs GTest into build directory
if(GTest_ARGUMENTS)
    cmake_parse_arguments(GTest
        "GIT_REPOSITORY;TIMEOUT"
        ""
        ${GTest_ARGUMENTS}
    )
endif()

set(arguments GIT_REPOSITORY)
if(GTest_GIT_REPOSITORY)
    list(APPEND arguments ${GTest_GIT_REPOSITORY})
else()
    list(APPEND arguments https://github.com/google/googletest)
endif()
if(GTest_TIMEOUT)
    list(APPEND arguments TIMEOUT ${GTest_TIMEOUT})
else()
    list(APPEND arguments TIMEOUT 10)
endif()

# write subset of variables to cache for gtest to use
include(PassonVariables)
passon_variables(GTest
    FILENAME "${EXTERNAL_ROOT}/src/GTestVariables.cmake"
    PATTERNS
        "CMAKE_[^_]*_R?PATH"
        "CMAKE_C_.*"
        "CMAKE_CXX_.*"
)

set(cmake_args -DBUILD_SHARED_LIBS=OFF -Dgtest_force_shared_crt=ON)
if(MINGW)
  list(APPEND cmake_args -Dgtest_disable_pthreads=ON)
endif()

ExternalProject_Add(
    Lookup-GTest
    PREFIX "${EXTERNAL_ROOT}"
    ${arguments}
    # Force separate output paths for debug and release builds to allow easy
    # identification of correct lib in subsequent TARGET_LINK_LIBRARIES commands
    CMAKE_ARGS
        -C "${EXTERNAL_ROOT}/src/GTestVariables.cmake"
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT}
        ${cmake_args}

    # Wrap download, configure and build steps in a script to log output
    UPDATE_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

add_recursive_cmake_step(Lookup-GTest DEPENDEES install)
