# Installs gabime/spdlog into build directory
#
# - GIT_REPOSITORY: defaults to origin gabime repo on github
# - GIT_TAG: defaults to master

if(spdlog_ARGUMENTS)
  cmake_parse_arguments(spdlog "" "GIT_REPOSITORY;GIT_TAG" ""
    ${spdlog_ARGUMENTS})
endif()
if(NOT spdlog_GIT_REPOSITORY)
  set(spdlog_GIT_REPOSITORY https://github.com/gabime/spdlog.git)
endif()
if(NOT spdlog_GIT_TAG)
  set(spdlog_GIT_TAG master)
endif()
ExternalProject_Add(
    Lookup-spdlog
    GIT_REPOSITORY ${spdlog_GIT_REPOSITORY}
    GIT_TAG ${spdlog_GIT_TAG}
    PREFIX "${EXTERNAL_ROOT}"
    # Force separate output paths for debug and release builds to allow easy
    # identification of correct lib in subsequent TARGET_LINK_LIBRARIES commands
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND
      ${CMAKE_COMMAND} -E copy_directory
          ${EXTERNAL_ROOT}/src/Lookup-spdlog/include
          ${EXTERNAL_ROOT}/include

    # Wrap download, configure and build steps in a script to log output
    UPDATE_COMMAND ""
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
)

add_recursive_cmake_step(Lookup-spdlog DEPENDEES install)
