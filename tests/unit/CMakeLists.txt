include_directories(${CMAKE_BINARY_DIR}/include)

add_definitions(-DBOOST_TEST_DYN_LINK)
#find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/include)

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

# Meshes
file(GLOB_RECURSE TEST_MESHES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    *.msh)
if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    foreach (f ${TEST_MESHES})
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${f}
                       ${CMAKE_CURRENT_BINARY_DIR}/${f}
                       COPYONLY)
    endforeach ()
endif ()
