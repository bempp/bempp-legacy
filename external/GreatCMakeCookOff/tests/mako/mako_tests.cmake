find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()

find_package(CoherentPython)

# Check location of binaries in build
set(origin "${PROJECT_BINARY_DIR}/../mako_build")
set(paths_exist "__init__.py" "other.py")
foreach(filename ${paths_exist})
    set(pathname "${origin}/python_package/makoed/${filename}")
    if(NOT EXISTS "${pathname}")
        message(FATAL_ERROR "Path ${pathname} not found")
    endif()
    set(pathname "${origin}/install/${filename}")
    if(NOT EXISTS "${pathname}")
        message(FATAL_ERROR "Path ${pathname} not found")
    endif()
endforeach()

# Check location does not exist
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} __init__.py
    WORKING_DIRECTORY "${origin}/python_package/makoed"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error
)
if(NOT result EQUAL 0)
    message("output: ${output}\n")
    message("error: ${error}\n")
    message("result: ${result}\n")
    message(FATAL_ERROR "Could not run mako")
endif()
