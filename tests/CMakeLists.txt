add_definitions(-DBOOST_TEST_DYN_LINK)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/lib ${Boost_INCLUDE_DIRS})
link_directories(${PROJECT_SOURCE_DIR}/lib ${Boost_LIBRARY_DIRS})

# put the names of all *.cpp files in variable test_source
file(GLOB_RECURSE test_sources *.cpp)
add_executable(run_tests
  ${test_sources}
  )
target_link_libraries(run_tests bempp ${Boost_LIBRARIES} dunecommon dunegrid)

add_custom_target(test
  ${CMAKE_CURRENT_BINARY_DIR}/run_tests
  COMMENT "Run unit tests" VERBATIM
  )
