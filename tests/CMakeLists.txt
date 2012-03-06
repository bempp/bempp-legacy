add_definitions(-DBOOST_TEST_DYN_LINK)
#find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/lib ${Boost_INCLUDE_DIR})
find_library(BOOST_UNIT_TEST 
	    NAMES libboost_unit_test_framework.so libboost_unit_test_framework-mt.so
            HINTS ${BOOST_LIBRARY_DIR} REQUIRED)

file(GLOB_RECURSE TEST_SOURCES *.cpp)
file(GLOB_RECURSE TEST_HEADERS *.hpp)
add_executable(run_tests
  ${TEST_SOURCES} ${TEST_HEADERS}
  )
target_link_libraries(run_tests bempp ${BOOST_UNIT_TEST})

add_custom_target(test
  ${CMAKE_CURRENT_BINARY_DIR}/run_tests
  COMMENT "Run unit tests" VERBATIM
  )
