# Immediately creates a directory
function(mkdir directory)
    if(NOT EXISTS "${directory}")
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E make_directory ${directory}
            OUTPUT_QUIET
        )
    endif()
endfunction()

# Immediately creates a symbolic link between two files
function(symlink FROM TO)
    if(NOT EXISTS "${FROM}")
        return()
    endif()
    if(EXISTS "${TO}")
        return()
    endif()
    if(WIN32)
        set(linkme mklink)
        if(IS_DIRECTORY "${FROM}")
          set(argument "/d")
        else()
          set(argument "")
        endif()
    else()
        set(linkme "ln")
        set(argument "-s")
    endif()
    get_filename_component(WD "${TO}" PATH)
    get_filename_component(TO "${TO}" NAME)
    execute_process(COMMAND ${linkme} ${argument} ${FROM} ${TO}
        WORKING_DIRECTORY ${WD}
    )
endfunction()

# Adds to a variety of environment variables
function(add_to_envvar VARIABLE PATH)
    include(CMakeParseArguments)
    cmake_parse_arguments(envvar "PREPEND" "OS" "" ${ARGN})
    if(envvar_OS AND NOT ${${envvar_OS}})
        return()
    endif()
    if("$ENV{${VARIABLE}}" STREQUAL "")
        set(separator "")
    elseif(WIN32)
        set(separator ";")
    else()
        set(separator ":")
    endif()
    if(envvar_PREPEND)
        set(ENV{${VARIABLE}} "${PATH}${separator}$ENV{${VARIABLE}}")
    else()
        set(ENV{${VARIABLE}} "$ENV{${VARIABLE}}${separator}${PATH}")
    endif()
endfunction()


