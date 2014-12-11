# Collects all the directories with *.dylibs and *.so in install manifest
# Then adds rpaths from all to all
# And changes tbb libraries
# This brute force approach ensures that rpath is set for external libraries as well.
unset(shared_objects)
unset(directories)
foreach(object ${CMAKE_INSTALL_MANIFEST_FILES})
    if(NOT IS_SYMLINK "${object}" AND NOT IS_DIRECTORY "${object}")
        if("${object}" MATCHES ".*\\.dylib$" OR "${object}" MATCHES ".*\\.so\\..*$")
            list(APPEND shared_objects "${object}")
            get_filename_component(directory "${object}" PATH)
            list(APPEND directories "${directory}")
        elseif("${object}" MATCHES ".*\\.so$")
            list(APPEND shared_objects "${object}")
            get_filename_component(filename "${object}" NAME)
            if(NOT "${filename}" MATCHES "_.*")
                get_filename_component(directory "${object}" PATH)
                list(APPEND directories "${directory}")
            endif()
        endif()
    endif()
endforeach()
list(REMOVE_DUPLICATES shared_objects)
list(REMOVE_DUPLICATES directories)

# Gets list of libraries that need @rpath added to on Darwin
unset(rpathless_libraries)
if("@CMAKE_SYSTEM_NAME@" STREQUAL "Darwin")
    set(libraries
        "@TBB_LIBRARY@" "@TBB_LIBRARY_DEBUG@"
        "@TBB_MALLOC_LIBRARY@" "@TBB_MALLOC_LIBRARY_DEBUG@"
    )
    set(NumpyMKL_FOUND @NumpyMKL_FOUND@)
    if(NumpyMKL_FOUND)
        list(APPEND libraries @NumpyMKL_LIBRARIES@)
    endif()
    foreach(library ${libraries})
        if(library AND EXISTS "${library}")
            get_filename_component(filename "${library}" NAME)
            list(APPEND rpathless_libraries "${filename}")
        endif()
    endforeach()
endif()

# Set up a macro to change rpath on given machine
if("@CMAKE_SYSTEM_NAME@" STREQUAL "Darwin")
    function(change_object_rpath object)
        foreach(directory ${ARGN})
            execute_process(COMMAND
                install_name_tool -add_rpath ${directory} ${object}
                ERROR_QUIET
            )
        endforeach()
    endfunction()
elseif("@CMAKE_SYSTEM_NAME@" STREQUAL "Linux")
    function(change_object_rpath object)
        get_filename_component(object_name "${object}" NAME_WE)
        if("${object_name}" MATCHES "libtbb.*")
            return() # libtbbs have not elf section, apparently.
        endif()
        get_filename_component(object_directory "${object}" PATH)
        set(rpath_string "$ORIGIN")
        foreach(directory ${ARGN})
            file(RELATIVE_PATH relative "${object_directory}" "${directory}")
            set(rpath_string "${rpath_string}:$ORIGIN/${relative}")
        endforeach()
        message(STATUS "RPath patching ${object}")
        execute_process(COMMAND
            @PATCHELF_EXECUTABLE@ --set-rpath ${rpath_string} ${object}
            RESULT_VARIABLE patchelf_result
        )
    endfunction()
else()
    message(FATAL_ERROR "Don't know how to change RPATH on @CMAKE_SYSTEM_NAME@")
endif()

# Loop over each shared object and add rpaths and change tbbs
foreach(library ${shared_objects})
    get_filename_component(current "${library}" PATH)
    set(current_directories ${directories})
    list(REMOVE_ITEM current_directories "${current}")
    change_object_rpath("${library}" ${current_directories})

    foreach(rpathless ${rpathless_libraries})
        execute_process(COMMAND
            install_name_tool -change ${rpathless} @rpath/${rpathless} ${library}
            ERROR_QUIET
        )
    endforeach()
endforeach()
