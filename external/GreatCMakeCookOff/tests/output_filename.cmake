find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()
include(ConfigureFiles)


function(check_same actual expected)
    get_filename_component(actual "${actual}" ABSOLUTE)
    get_filename_component(expected "${expected}" ABSOLUTE)
    if(NOT "${expected}" STREQUAL "${actual}")
        message(FATAL_ERROR "Unexpected output file: ${actual} vs ${expected}")
    endif()
endfunction()

output_filename("hello" output)
check_same("${output}" "${CMAKE_CURRENT_BINARY_DIR}/hello")

output_filename("hello/world" output)
check_same("${output}" "${CMAKE_CURRENT_BINARY_DIR}/hello/world")

output_filename("hello/world" output "${CMAKE_CURRENT_SOURCE_DIR}/there")
check_same("${output}" "${CMAKE_CURRENT_SOURCE_DIR}/there/hello/world")
