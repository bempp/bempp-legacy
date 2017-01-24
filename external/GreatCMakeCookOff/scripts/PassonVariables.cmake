# A script to write cache variables to file
# Makes it easy to include a subset of cached variables in external projects
function(passon_variables PACKAGE)
    include(CMakeParseArguments)
    cmake_parse_arguments(passon
        "APPEND;PUBLIC" "FILENAME" "PATTERNS;ALSOADD"
        ${ARGN}
    )
    if(NOT passon_FILENAME AND EXTERNAL_ROOT)
        set(passon_FILENAME "${EXTERNAL_ROOT}/src/${PACKAGE}.cmake")
    elseif()
        set(passon_FILENAME "${CMAKE_BINARY_DIR}/CMakeFiles/${PACKAGE}.cmake")
    endif()
    get_cmake_property(all_cached_variables CACHE_VARIABLES)
    if(passon_PUBLIC)
        set(oldvarlist ${all_cached_variables})
        set(all_cached_variables "")
        foreach(variable ${oldvarlist})
            get_property(type CACHE ${variable} PROPERTY TYPE)
            if(NOT "${type}" STREQUAL "INTERNAL")
                list(APPEND all_cached_variables ${variable})
            endif()
        endforeach()
    endif()
    set(setters "")
    foreach(variable ${all_cached_variables})
        set(does_match False)
        foreach(pattern ${passon_PATTERNS})
            if(variable MATCHES "${pattern}")
                set(does_match True)
                break()
            endif()
        endforeach()
        if(does_match)
            get_property(type CACHE ${variable} PROPERTY TYPE)
            get_property(help CACHE ${variable} PROPERTY HELPSTRING)
            set(setters
              "${setters}\nset(${variable} \"${${variable}}\" CACHE ${type} \"${help}\")"
            )
        endif()
    endforeach()
    if(NOT passon_APPEND)
      file(WRITE "${passon_FILENAME}" "# pre-cached variables for ${PACKAGE}")
    endif()
    file(APPEND "${passon_FILENAME}"
        "${setters}\n"
        "\n# Explicitely added lines\n"
        ${passon_ALSOADD}
        "\n# End of passon_variables\n"
    )
endfunction()
