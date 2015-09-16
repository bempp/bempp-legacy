include(ConfigureFiles)

function(copy_files_on_change)
    # Parses arguments
    cmake_parse_arguments(_copy_files
        ""
        ";OUTPUT_FILES;DESTINATION;TARGETNAME"
        "CMDLINE;DEPENDENCIES;GLOB;DEPENDS"
        ${ARGN}
    )
    unset(sources)
    if(NOT "${_copy_files_GLOB}" STREQUAL "")
        file(GLOB sources ${_copy_files_GLOB})
    endif()
    list(APPEND sources ${_copy_files_UNPARSED_ARGUMENTS})
    list(REMOVE_DUPLICATES sources)
    if("${sources}" STREQUAL "")
        return()
    endif()
    set(destination "${_copy_files_DESTINATION}")
    if("${destination}" STREQUAL "")
        set(destination "${CMAKE_CURRENT_BINARY_DIR}")
    else()
        file(MAKE_DIRECTORY "${destination}")
    endif()

    unset(copied_files)
    foreach(filename ${sources})
        output_filename("${filename}" output "${destination}")
        get_filename_component(abspath "${filename}" ABSOLUTE)
        add_custom_command(
            OUTPUT ${output}
            COMMAND ${CMAKE_COMMAND} -E copy ${abspath} ${output}
            MAIN_DEPENDENCY ${abspath}
            COMMENT "Copying file ${filename}")
        list(APPEND copied_files "${output}")
    endforeach()

    set(${_copy_files_OUTPUT_FILES} ${copied_files} PARENT_SCOPE)
    # Add custom target. Makes it easier to handle dependencies.
    if(NOT "${_copy_files_TARGETNAME}" STREQUAL "")
        add_custom_target(${_copy_files_TARGETNAME} DEPENDS ${copied_files})
    endif()
endfunction()
