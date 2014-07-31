if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    include(TargetCopyFiles)
    add_custom_target(copy_python_integration_tests ALL)
    add_copy_files(copy_python_integration_tests
        GLOB  *.py
        DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
    )
endif ()
