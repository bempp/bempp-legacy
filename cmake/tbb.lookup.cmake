# Not clear that parsing works in find_package. So setting variable directly.
find_package(TBB)

if(NOT TBB_FOUND)
  message(STATUS "Threading Building Blocks not found. Will attempt to download it.")
  
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(
      TBB_URL
      http://threadingbuildingblocks.org/sites/default/files/software_releases/mac/tbb42_20131003oss_osx.tgz
    )
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(
      TBB_URL
      http://threadingbuildingblocks.org/sites/default/files/software_releases/linux/tbb42_20131003oss_lin.tgz
    )
  else()
    message(FATAL_ERROR "Automatic install of TBB in windows not implemented")
  endif()
  ExternalProject_Add(
      tbb
      PREFIX ${EXTERNAL_ROOT}
      # Downloads boost from url -- much faster than svn
      URL ${TBB_URL}
      BUILD_IN_SOURCE 1
      PATCH_COMMAND
       ${CMAKE_COMMAND} -E copy_if_different
                        ${CMAKE_SOURCE_DIR}/cmake/tbb.install.cmake
                        ${EXTERNAL_ROOT}/src/tbb/CMakeLists.txt
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT}
                 -DROOT=${PROJECT_SOURCE_DIR}
      LOG_DOWNLOAD ON
      LOG_CONFIGURE ON
      LOG_BUILD ON
  )
  # Rerun cmake to capture new boost install
  add_recursive_cmake_step(tbb DEPENDEES install)
    
endif()
