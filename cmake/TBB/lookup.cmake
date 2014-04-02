# Not clear that parsing works in find_package. So setting variable directly.
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(
    TBB_URL
    http://threadingbuildingblocks.org/sites/default/files/software_releases/mac/tbb42_20131003oss_osx.tgz
  )
  set(TBB_URL_MD5 bb714437fc5443e196309299ea473570)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(
    TBB_URL
    http://threadingbuildingblocks.org/sites/default/files/software_releases/linux/tbb42_20131003oss_lin.tgz
  )
  set(TBB_URL_MD5 8e43d060fb0303f2437b8e84c7762582)
else()
  message(FATAL_ERROR "Automatic install of TBB in windows not implemented")
endif()
ExternalProject_Add(
    TBB
    PREFIX ${EXTERNAL_ROOT}
    # Downloads boost from url -- much faster than svn
    URL ${TBB_URL}
    URL_HASH MD5=${TBB_URL_MD5}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND
     ${CMAKE_COMMAND} -E copy_if_different
                      ${CURRENT_LOOKUP_DIRECTORY}/install.cmake
                      ${EXTERNAL_ROOT}/src/tbb/CMakeLists.txt
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT}
               -DROOT=${PROJECT_SOURCE_DIR}
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)
# Rerun cmake to capture new TBB install
add_recursive_cmake_step(TBB TBB_FOUND DEPENDEES install)
