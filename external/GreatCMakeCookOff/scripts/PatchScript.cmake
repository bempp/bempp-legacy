# Creates a script that applies patches
# This script does not fail on errors. It ensures that we can "apply" patches over and over again.
include(CMakeParseArguments)
find_program(PATCH_EXECUTABLE patch)
find_program(BASH_EXECUTABLE bash)
function(create_patch_script NAME OUTVAR)
    if(NOT PATCH_EXECUTABLE)
        message(FATAL_ERROR "Could not find the patch program")
    endif()
    cmake_parse_arguments(patcher
        ""
        "CMDLINE;WORKING_DIRECTORY"
        ""
        ${ARGN}
    )
    if(NOT patcher_CMDLINE)
        set(patcher_CMDLINE "")
    endif()
    if(NOT patcher_WORKING_DIRECTORY)
        set(patcher_WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    endif()

    # Create patcher script
    set(script_file "${PROJECT_BINARY_DIR}/CMakeFiles/patches/noperms/${NAME}.sh")
    file(WRITE "${script_file}"
        "#!${BASH_EXECUTABLE}\n"
        "cd ${patcher_WORKING_DIRECTORY}\n"
    )
    foreach(filename ${patcher_UNPARSED_ARGUMENTS})
        get_filename_component(filename "${filename}" ABSOLUTE)
        file(APPEND "${script_file}"
            "${PATCH_EXECUTABLE} -N ${patcher_CMDLINE} < ${filename}\n"
        )
    endforeach()
    file(APPEND "${script_file}" "true\n")

    file(COPY "${script_file}"
        DESTINATION "${PROJECT_BINARY_DIR}/CMakeFiles/patches/"
        FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
    )

    set(${OUTVAR} "${PROJECT_BINARY_DIR}/CMakeFiles/patches/${NAME}.sh" PARENT_SCOPE)
endfunction()
