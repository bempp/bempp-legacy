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
  
  # Downloads boost from svn
  ExternalProject_Add(
      boost
      PREFIX ${EXTERNAL_ROOT}
      SVN_REPOSITORY http://svn.boost.org/svn/boost/tags/release/Boost_1_55_0
      TIMEOUT 10
      UPDATE_COMMAND ""
      CONFIGURE_COMMAND ""
      BUILD_COMMAND  ""
      INSTALL_COMMAND ""
      # Wrap download, configure and build steps in a script to log output
      LOG_DOWNLOAD ON
      LOG_CONFIGURE ON
      LOG_BUILD ON
  )
  # create bjam
  ExternalProject_Add_Step(
    boost Configure
    DEPENDEES download
    COMMAND ./bootstrap.sh
    WORKING_DIRECTORY ${EXTERNAL_ROOT}/src/boost
  )
  # build
  ExternalProject_Add_Step(
    boost Build
    DEPENDEES Configure
    COMMAND ./b2 link=static variant=release --with-test 
    WORKING_DIRECTORY ${EXTERNAL_ROOT}/src/boost
  )
  # install
  ExternalProject_Add_Step(
    boost Install
    DEPENDEES Build
    COMMAND ./b2 link=static variant=release --with-test --prefix=${EXTERNAL_ROOT} install
    WORKING_DIRECTORY ${EXTERNAL_ROOT}/src/boost
  )
  # Rerun cmake to capture new boost install
  add_rerun_cmake_step(boost DEPENDEES Install)
    
endif()

message(STATUS "[Boost unit-test] ${BOOST_UNIT_TEST_LIB}")
