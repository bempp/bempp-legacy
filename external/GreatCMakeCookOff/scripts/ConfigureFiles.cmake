include(CMakeParseArguments)
function(output_filename filename OUTPUT)
    list(LENGTH ARGN remaining_args)
    if(remaining_args EQUAL 0)
        set(destination "${CMAKE_CURRENT_BINARY_DIR}")
    elseif(remaining_args EQUAL 1)
        set(destination "${ARGN}")
    else()
        message(FATAL_ERROR "Too many arguments to output_filename: ${ARGN}")
    endif()

    get_filename_component(filename "${filename}" ABSOLUTE)
    file(RELATIVE_PATH relfile "${CMAKE_CURRENT_SOURCE_DIR}" "${filename}")

    if("${relfile}" MATCHES "\\.\\./")
        file(RELATIVE_PATH relfile "${destination}" "${filename}")
        if("${relfile}" MATCHES "\\.\\./")
            message(FATAL_ERROR "File ${filename} is not in "
                "destination or source directory or subdirectory.")
        endif()
    endif()
    set(${OUTPUT} "${destination}/${relfile}" PARENT_SCOPE)
endfunction()

macro(configure_files)
    # Parses arguments
    cmake_parse_arguments(_cf
        "" "OUTPUT_FILES;DESTINATION" "GLOB"
        ${ARGN}
    )

    unset(_cf_sources)
    if(NOT "${_mako_GLOB}" STREQUAL "")
        file(GLOB _cf_sources ${_cf_GLOB})
    endif()
    list(APPEND _cf_sources ${_cf_UNPARSED_ARGUMENTS})
    list(REMOVE_DUPLICATES _cf_sources)
    if("${_cf_sources}" STREQUAL "")
        return()
    endif()

    unset(_cf_configured_files)
    foreach(filename ${_cf_sources})
        output_filename("${filename}" output "${_cf_DESTINATION}")
        string(REGEX REPLACE "(.*)\\.in(\\..*)" "\\1\\2" output "${output}")

        configure_file("${filename}" "${output}" @ONLY)
        list(APPEND _cf_configured_files "${output}")
    endforeach()

    if(NOT _cf_OUTPUT_FILES STREQUAL "")
        set(${_cf_OUTPUT_FILES} ${_cf_configured_files})
    endif()
endmacro()

