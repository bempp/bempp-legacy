# Scripts that modify LD_LIBRARY_PATH and such

# Adds a single path to a path-file
function(_add_to_a_path_single THISFILE path)
    if("${path}" STREQUAL "" OR "${path}" MATCHES "NOTFOUND")
        return()
    endif()
    # If quacks like a library, get directory where it resides
    get_filename_component(extension "${path}" EXT)
    if("${extension}" MATCHES "\\.so.*" OR "${extension}" MATCHES "\\.dylib")
        get_filename_component(path "${path}" PATH)
    elseif("${extension}" MATCHES "\\.a")
        return() # Archive are not dynamic, no need to add to rpath.
    endif()
    # Makes it an absolute path
    get_filename_component(path "${path}" ABSOLUTE)
    # Add to path file if not there yet
    if(NOT EXISTS "${THISFILE}")
        file(WRITE "${THISFILE}" "${path}\n")
    else()
        file(STRINGS "${THISFILE}" ALLPATHS)
        list(FIND ALLPATHS "${path}" INDEX)
        if(INDEX EQUAL -1)
            file(APPEND "${THISFILE}" "${path}\n")
        endif()
    endif()
endfunction()
# Adds many paths to a path file
function(_add_to_a_path THISFILE)
    foreach(path ${ARGN})
        _add_to_a_path_single("${THISFILE}" "${path}")
    endforeach()
endfunction()

function(add_to_ld_path)
    unset(ldpaths)
    foreach(directory ${ARGN})
        get_filename_component(directory "${directory}" ABSOLUTE)
        if(NOT directory MATCHES "^\\/System\\/")
            list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
                "${directory}" is_system_dir)
            if("${is_system_dir}" STREQUAL "-1")
                list(APPEND ldpaths "${directory}")
            endif()
        endif()
    endforeach()
    _add_to_a_path("${PROJECT_BINARY_DIR}/paths/ldpaths" ${ldpaths})
endfunction()
function(add_to_python_path)
    _add_to_a_path("${PROJECT_BINARY_DIR}/paths/pypaths.pth" ${ARGN})
endfunction()

# Gets list of eggs in given directories
function(_get_python_eggs OUTVAR)
    unset(patterns)
    foreach(directory ${ARGN})
        if(NOT "${directory}" MATCHES ".*egg")
            list(APPEND patterns "${directory}/*.egg")
        else()
            list(APPEND patterns "${directory}")
        endif()
    endforeach()
    if(patterns)
        file(GLOB directories ${patterns})
        set(${OUTVAR} ${directories} PARENT_SCOPE)
    endif()
endfunction()

# Add eggy directories to python path
function(add_python_eggs)
    cmake_parse_arguments(apeg "" "" "EXCLUDE;INCLUDE" ${ARGN})
    _get_python_eggs(included_eggs ${apeg_UNPARSED_ARGUMENTS} ${apeg_INCLUDE})
    _get_python_eggs(excluded_eggs ${apeg_EXCLUDE})
    if(excluded_eggs)
        list(REMOVE_ITEM included_eggs ${excluded_eggs})
    endif()
    add_to_python_path(${included_eggs})
endfunction()

if(NOT UNIX)
    function(create_environment_script caller location)
        message(FATAL_ERROR "Environment scripts not implemented "
            "on non-UNIX systems")
    endfunction()
    return()
endif()

get_filename_component(_PATH_TO_LOCALBASH_IN
    "${CMAKE_CURRENT_LIST_DIR}/localbash.in.sh"
    ABSOLUTE
)
find_program(BASH_EXECUTABLE bash)
find_program(ENV_EXECUTABLE env)
include(CMakeParseArguments)

function(create_environment_script)
    cmake_parse_arguments(env 
        "PYTHON"
        "SCRIPT;PATH;EXECUTABLE;WORKING_DIRECTORY"
        "" ${ARGN}
    )
    if(NOT env_PATH)
        set(env_PATH "${CMAKE_CURRENT_BINARY_DIR}/envscript.sh")
    endif()
    if(NOT env_EXECUTABLE)
        set(env_EXECUTABLE "")
    endif()
    # used in the configured script: if set, modifies python path
    if(NOT env_PYTHON)
        set(env_PYTHON "")
    endif()
    if(NOT env_SCRIPT)
        set(env_SCRIPT "${_PATH_TO_LOCALBASH_IN}")
    endif()

    get_filename_component(filename "${env_PATH}" NAME)
    get_filename_component(directory "${env_PATH}" PATH)
    configure_file("${env_SCRIPT}"
        "${PROJECT_BINARY_DIR}/CMakeFiles/${filename}"
        @ONLY
    )
    file(COPY "${PROJECT_BINARY_DIR}/CMakeFiles/${filename}"
        DESTINATION "${directory}"
        FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
    )
endfunction()
