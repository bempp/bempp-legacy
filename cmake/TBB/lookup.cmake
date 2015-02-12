# Not clear that parsing works in find_package. So setting variable directly.
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(
    TBB_URL
    https://www.threadingbuildingblocks.org/sites/default/files/software_releases/mac/tbb43_20150209oss_osx.tgz
  )
  set(TBB_URL_MD5 3e683c19792582b61382e0d760ea5db2)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(
    TBB_URL
    https://www.threadingbuildingblocks.org/sites/default/files/software_releases/linux/tbb43_20150209oss_lin.tgz
  )
  set(TBB_URL_MD5 8a13558019c9acfbab269e2b34f7fd27)
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
    URL_MD5 ${TBB_URL_MD5}
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
