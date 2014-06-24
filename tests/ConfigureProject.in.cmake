set(build_dir "@CMAKE_CURRENT_BINARY_DIR@/ExampleProject/build")
if(EXISTS "${build_dir}")
    message(STATUS "Removing previous build directory")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${build_dir}
        RESULT_VARIABLE result
    )
    if(NOT result STREQUAL 0)
        message(FATAL_ERROR "Could not remove directory")
    endif()
endif()
file(MAKE_DIRECTORY "${build_dir}")

# Make sure that some test project picks up some things we know from this
# build:
# - compiler should be the same
# - Bempp_DIR should point build dir, so that this build is picked up, rather
# than any other versions of bem++ that may be lying around.
file(WRITE "${build_dir}/CacheVar.cmake"
    "set(CMAKE_CXX_COMPILER \"@CMAKE_CXX_COMPILER@\" CACHE PATH \"\")\n"
    "set(Bempp_DIR \"@PROJECT_BINARY_DIR@\" CACHE PATH \"\")\n"
)

execute_process(
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_INCLUSION_UNIT_TESTS=ON
        -CCacheVar.cmake ..
    WORKING_DIRECTORY "${build_dir}"
    RESULT_VARIABLE result
)

if(NOT result STREQUAL 0)
    message("error: ${error}\n")
    message("output: ${output}\n")
    message(FATAL_ERROR "Could not configure test project in ${build_dir}")
endif()
