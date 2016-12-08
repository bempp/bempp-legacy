include(ConfigureFiles)

function(mako_files)
    # Parses arguments
    cmake_parse_arguments(_mako
        ""
        "MAKO_SCRIPT;OUTPUT_FILES;DESTINATION;TARGETNAME"
        "CMDLINE;DEPENDENCIES;GLOB;DEPENDS"
        ${ARGN}
    )
    unset(sources)
    if(NOT "${_mako_GLOB}" STREQUAL "")
        file(GLOB sources ${_mako_GLOB})
    endif()
    list(APPEND sources ${_mako_UNPARSED_ARGUMENTS})
    list(REMOVE_DUPLICATES sources)
    if("${sources}" STREQUAL "")
        return()
    endif()
    set(destination "${_mako_DESTINATION}")
    if("${destination}" STREQUAL "")
        set(destination "${CMAKE_CURRENT_BINARY_DIR}")
    else()
        file(MAKE_DIRECTORY "${destination}")
    endif()

    set(local_python "${LOCAL_PYTHON_EXECUTABLE}")
    if("${local_python}" STREQUAL "")
        set(local_python ${PYTHON_EXECUTABLE})
    endif()
    set(mako_script "${_mako_MAKO_SCRIPT}")
    if("${mako_script}" STREQUAL "")
        if("${mako_SCRIPT}" STREQUAL "")
            message(FATAL_ERROR "Mako render script not defined")
        endif()
        set(mako_script "${mako_SCRIPT}")
    endif()

    unset(makoed_files)
    foreach(filename ${sources})
        output_filename("${filename}" output "${destination}")
        string(REGEX REPLACE "(.*)\\.mako(\\..*)" "\\1\\2" output "${output}")
        get_filename_component(abspath "${filename}" ABSOLUTE)

        add_custom_command(
            OUTPUT ${output}
            COMMAND
                ${local_python} -B ${mako_script}
                    ${abspath} ${_mako_CMDLINE} > ${output}
            DEPENDS "${abspath}" "${mako_script}" ${_mako_DEPENDS}
            COMMENT "Mako-ing file ${filename}"
        )
        list(APPEND makoed_files "${output}")
    endforeach()

    set(${_mako_OUTPUT_FILES} ${makoed_files} PARENT_SCOPE)
    # Add custom target. Makes it easier to handle dependencies.
    if(NOT "${_mako_TARGETNAME}" STREQUAL "")
        add_custom_target(${_mako_TARGETNAME} DEPENDS ${makoed_files})
    endif()
endfunction()

