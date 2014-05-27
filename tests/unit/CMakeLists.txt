include(BoostUnitTestFlags)
if(Boost_UNIT_TEST_CXXFLAGS)
    add_definitions(${Boost_UNIT_TEST_CXXFLAGS})
endif()

file(GLOB_RECURSE TEST_SOURCES *.cpp)
file(GLOB_RECURSE TEST_HEADERS *.hpp)
add_executable(run_tests ${TEST_SOURCES} ${TEST_HEADERS})
depends_on_lookups(run_tests)

target_link_libraries(run_tests bempp)
if(Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG)
    target_link_libraries(run_tests
        debug ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG}
    )
endif()
if(Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE)
    target_link_libraries(run_tests
        optimized ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE}
    )
endif()

add_test(NAME unittests
    COMMAND ${PROJECT_BINARY_DIR}/run_cpptests.sh $<TARGET_FILE:run_tests>
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
if(NOT "$ENV{LD_LIBRARY_PATH}" STREQUAL "")
    list(APPEND ENVIRONMENT_VARIABLES
        "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}")
endif()
if(NOT "$ENV{DYLD_LIBRARY_PATH}" STREQUAL "")
    list(APPEND ENVIRONMENT_VARIABLES
        "DYLD_LIBRARY_PATH=$ENV{DYLD_LIBRARY_PATH}")
endif()
set_tests_properties(unittests
    PROPERTIES ENVIRONMENT "${ENVIRONMENT_VARIABLES}"
    LABELS "cpp;unit"
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
