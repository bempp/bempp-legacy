# Not clear that parsing works in find_package. So setting variable directly.
set(Boost_FIND_COMPONENTS unit_test_framework)
find_package(boost COMPONENTS ${Boost_FIND_COMPONENTS})

if(Boost_FOUND)
  if(CMAKE_BUILD_TYPE STREQUAL "Release" OR MAKE_BUILD_TYPE STREQUAL "RelWithDebInfo") 
    set(BOOST_UNIT_TEST_LIB ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})
  else()
    set(BOOST_UNIT_TEST_LIB ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG})
  endif()

  set(
    BOOST_UNIT_TEST_LIB ${BOOST_UNIT_TEST_LIB} 
    CACHE INTERNAL
    "Path to unit test framework"
  )
  set(
    BOOST_INCLUDE_DIR ${Boost_INCLUDE_DIR}
    CACHE INTERNAL
    "Path to boost include directory"
  )
else()
  message(STATUS "Boost not found. Will attempt to download it.")
  
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
  add_recursive_cmake_step(Boost Boost_FOUND DEPENDEES install)
    
endif()
