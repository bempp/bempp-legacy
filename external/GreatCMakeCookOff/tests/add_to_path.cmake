find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()
include(EnvironmentScript)

## Delete stuff from previous test
foreach(filename ldpaths pypaths.pth)
    set(filename "${PROJECT_BINARY_DIR}/paths/${filename}")
    if(EXISTS "${filename}")
        execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${filename})
    endif()
endforeach()

## First checks ldpaths: adds a relative and an absolute path
add_to_ld_path("relative/path" "/absolute/path")
# Check file exists
if(NOT EXISTS "${PROJECT_BINARY_DIR}/paths/ldpaths")
    message(FATAL_ERROR "ldpaths file does not exist")
endif()

file(STRINGS "${PROJECT_BINARY_DIR}/paths/ldpaths" lines)
list(LENGTH lines length)
if(NOT length EQUAL 2)
    message(FATAL_ERROR "Incorrect list length")
endif()
list(GET lines 0 line)
get_filename_component(abspath "relative/path" ABSOLUTE)
if(NOT "${line}" STREQUAL "${abspath}")
    message(FATAL_ERROR "Expected path not found in ldpaths")
endif()
list(GET lines 1 line)
if(NOT "${line}" STREQUAL "/absolute/path")
    message(FATAL_ERROR "Expected path not found in ldpaths")
endif()

# Adds same paths: nothing should happen
add_to_ld_path("/absolute/path" "relative/path")
file(STRINGS "${PROJECT_BINARY_DIR}/paths/ldpaths" lines)
list(LENGTH lines length)
if(NOT length EQUAL 2)
    message(FATAL_ERROR "Incorrect list length")
endif()

# Add same path as library: nothing should happen
add_to_ld_path("/absolute/path/libsomething.so")
file(STRINGS "${PROJECT_BINARY_DIR}/paths/ldpaths" lines)
list(LENGTH lines length)
if(NOT length EQUAL 2)
    message(FATAL_ERROR "Incorrect list length")
endif()

# Add different path as archive: nothing  should happen
add_to_ld_path("/absolute/other/path/libsomething.a")
file(STRINGS "${PROJECT_BINARY_DIR}/paths/ldpaths" lines)
list(LENGTH lines length)
if(NOT length EQUAL 2)
    message(FATAL_ERROR "Incorrect list length")
endif()

# Add different path as library
add_to_ld_path("/absolute/other/path/libsomething.so")
file(STRINGS "${PROJECT_BINARY_DIR}/paths/ldpaths" lines)
if(NOT "${lines}" STREQUAL "${abspath};/absolute/path;/absolute/other/path")
    message(FATAL_ERROR "Incorrect ld paths ${lines}")
endif()

# Add different path as library
add_to_ld_path("/other/path/libsomething.dylib")
file(STRINGS "${PROJECT_BINARY_DIR}/paths/ldpaths" lines)
if(NOT "${lines}" STREQUAL "${abspath};/absolute/path;/absolute/other/path;/other/path")
    message(FATAL_ERROR "Incorrect ld paths ${lines}")
endif()

## Then checks pypaths
add_to_python_path("relative/path" "/absolute/path")
if(NOT EXISTS "${PROJECT_BINARY_DIR}/paths/pypaths.pth")
    message(FATAL_ERROR "ldpaths file does not exist")
endif()
file(STRINGS "${PROJECT_BINARY_DIR}/paths/pypaths.pth" lines)
if(NOT "${lines}" STREQUAL "${abspath};/absolute/path")
    message(FATAL_ERROR "Incorrect py paths ${lines}")
endif()
