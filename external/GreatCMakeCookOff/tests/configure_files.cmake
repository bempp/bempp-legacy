find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()
include(ConfigureFiles)

file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/first.in.h"
    "hello @PROJECT_BINARY_DIR@")
file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/second.py"
    "hello @PROJECT_BINARY_DIR@")

configure_files(
    OUTPUT_FILES output
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/there"
    first.in.h
    second.py
)

foreach(filename first.h second.py)
    set(filename "${CMAKE_CURRENT_BINARY_DIR}/there/${filename}")
    if(NOT EXISTS "${filename}")
        message("${filename} not found")
    endif()
    list(FIND output "${filename}" found)
        message("${filename} not in output ${output}")
    if(found LESS 0)
        message(FATAL_ERROR "${filename} not in output")
    endif()
endforeach()

list(LENGTH output i)
if(NOT i EQUAL 2)
    message(FATAL_ERROR "Incorrect number of output files")
endif()
