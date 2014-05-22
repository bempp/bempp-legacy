message(STATUS "Installing Trilinos")
execute_process(
    COMMAND ${CMAKE_COMMAND}
        -DPyTrilinos_INSTALL_DIR=@install_location@
        -DPyTrilinos_INSTALL_PREFIX=@install_location@
        -DCMAKE_INSTALL_PREFIX=@install_location .
    WORKING_DIRECTORY "@EXTERNAL_ROOT@/src/Trilinos-build/"
    RESULT_VARIABLE result
    ERROR_VARIABLE error
)
execute_process(
    COMMAND ${CMAKE_COMMAND} --build . --target install
    WORKING_DIRECTORY "@EXTERNAL_ROOT@/src/Trilinos-build/"
    RESULT_VARIABLE result
    ERROR_VARIABLE error
)
if(NOT ${result} EQUAL 0)
    message("error: ${error}")
    message("error code: ${result}")
    message(FATAL_ERROR "Could not install Trilinos to @install_location@")
endif()
execute_process(
    COMMAND ${CMAKE_COMMAND} -C ${EXTERNAL_ROOT}/src/TrilinosVariables.cmake .
    WORKING_DIRECTORY "@EXTERNAL_ROOT@/src/Trilinos-build/"
    RESULT_VARIABLE result
    ERROR_VARIABLE error
)
