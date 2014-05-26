add_subdirectory(unit)
if (WITH_INTEGRATION_TESTS)
   add_subdirectory(integration)
endif ()

# Tries and creates an example project
if(NOT NOEXPORT)
    configure_file("ExampleProject_CMakeLists.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/ExampleProject/CMakeLists.txt"
        @ONLY
    )
    configure_file("ProjectInclusionUnitTests.in.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/ExampleProject/ProjectInclusionUnitTests.cmake"
        @ONLY
    )
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/ExampleProject/build")
    add_test(NAME configure_test_project
        COMMAND ${CMAKE_COMMAND} .. -DPROJECT_INCLUSION_UNIT_TESTS=ON
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/ExampleProject/build"
    )
    add_test(NAME build_test_project
        COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/ExampleProject/build"
    )
    set_tests_properties(build_test_project PROPERTIES DEPENDS configure_test_project)
    set_tests_properties(build_test_project configure_test_project
        PROPERTIES LABELS "cmake")
endif()
