# Not clear that parsing works in find_package. So setting variable directly.
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(
    TBB_URL
    http://www.threadingbuildingblocks.org/sites/default/files/software_releases/mac/tbb42_20140601oss_osx.tgz
  )
  set(TBB_URL_SHA256 
      63078757f3f65abb3ef568e75568e721fa83efdecad4d6a060b010255ddbdee3)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(
    TBB_URL
    http://www.threadingbuildingblocks.org/sites/default/files/software_releases/linux/tbb42_20140601oss_lin.tgz
  )
  set(TBB_URL_SHA256
      95aab3f82f435aff8aa1f63d1e3608e3181fc66d2d52d0c429873af8c69af882)
else()
  message(FATAL_ERROR "Automatic install of TBB in windows not implemented")
endif()

# Writes a script of cached variables
file(WRITE "${EXTERNAL_ROOT}/src/TBBVariables.cmake"
    "\nset(CMAKE_INSTALL_PREFIX \"${EXTERNAL_ROOT}\" CACHE PATH \"\")\n"
    "\nset(CMAKE_MODULE_PATH \"${CMAKE_MODULE_PATH}\" CACHE STRING \"\")\n"
    "\nset(ROOT \"${PROJECT_SOURCE_DIR}\" CACHE STRING \"\")\n"
)

ExternalProject_Add(
    TBB
    PREFIX ${EXTERNAL_ROOT}
    URL ${TBB_URL}
    URL_HASH SHA256=${TBB_URL_SHA256}
    PATCH_COMMAND 
        ${CMAKE_COMMAND} -E copy_if_different
            ${CURRENT_LOOKUP_DIRECTORY}/install.cmake
            ${EXTERNAL_ROOT}/src/TBB/CMakeLists.txt
    CMAKE_ARGS  -C "${EXTERNAL_ROOT}/src/TBBVariables.cmake"
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)


# Rerun cmake to capture new TBB install
add_recursive_cmake_step(TBB DEPENDEES install)
