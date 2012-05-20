include_directories(${CMAKE_BINARY_DIR}/include/common)

add_definitions(-DBOOST_TEST_DYN_LINK)
#find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/lib ${Boost_INCLUDE_DIR})


file(GLOB_RECURSE TEST_SOURCES *.cpp)
file(GLOB_RECURSE TEST_HEADERS *.hpp)
add_executable(run_tests
  ${TEST_SOURCES} ${TEST_HEADERS}
  )
  
target_link_libraries(run_tests bempp ${BOOST_UNIT_TEST_LIB})

add_custom_target(test
  ${CMAKE_CURRENT_BINARY_DIR}/run_tests
  COMMENT "Run unit tests" VERBATIM
  )
