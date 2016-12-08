function(test_exists test_name)
    execute_process(
        COMMAND ${CMAKE_CTEST_COMMAND} -N -R ${test_name}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../pytest_build/"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error
    )
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Could not run ctest command")
    endif()
    string(REGEX REPLACE ".*Total Tests: ([0-9]+).*" "\\1"
        nb_tests ${output})
    if(NOT nb_tests EQUAL 1)
        message(FATAL_ERROR "Could not find test ${test_name}")
    endif()
endfunction()

function(check_test_passes test_name expected_ntests)
    execute_process(
        COMMAND ${CMAKE_CTEST_COMMAND} -V -R ${test_name}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../pytest_build/"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_QUIET
    )
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "test ${test_name} did not pass - ${result}")
    endif()
    if(NOT "${output}" MATCHES "collected ${expected_ntests} items")
        message(STATUS "o: ${output}")
        message(FATAL_ERROR
            "Expected ${expected_ntests} tests in ${test_name}."
        )
    endif()
endfunction()
function(check_test_fails test_name expected_ntests)
    execute_process(
        COMMAND ${CMAKE_CTEST_COMMAND} -V -R ${test_name}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../pytest_build/"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_QUIET
    )
    if(result EQUAL 0)
        message(FATAL_ERROR "test ${test_name} did not fail")
    endif()
    if(NOT "${output}" MATCHES "collected ${expected_ntests} items")
        message(STATUS "o: ${output}")
        message(FATAL_ERROR
            "Expected ${expected_ntests} in ${test_name}, got ${ntests}."
        )
    endif()
endfunction()

test_exists(hackage.this)
test_exists(hackage.that)
test_exists(hackage.cmdl)
test_exists(hackage.fails.cmdl)
test_exists(hackage.cython)
check_test_passes(hackage.this 2)
check_test_fails(hackage.that 1)
check_test_passes(hackage.cmdl 1)
check_test_fails(hackage.fails.cmdl 1)
check_test_passes(hackage.cython 2)
