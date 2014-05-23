# Collects all the directories with *.dylibs and *.so in install manifest
# Then adds rpaths from all to all
# And changes tbb libraries
# This brute force approach ensures that rpath is set for external libraries as well.
unset(shared_objects)
unset(directories)
foreach(object ${CMAKE_INSTALL_MANIFEST_FILES})
    if("${object}" MATCHES ".*\\.dylib$")
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
endforeach()
list(REMOVE_DUPLICATES shared_objects)
list(REMOVE_DUPLICATES directories)

# Gets list of tbbs on Darwin
unset(tbb_names)
if("@CMAKE_SYSTEM_NAME@" STREQUAL "Darwin")
    foreach(tbb "@TBB_LIBRARY@" "@TBB_LIBRARY_DEBUG@"
        "@TBB_MALLOC_LIBRARY@" "@TBB_MALLOC_LIBRARY_DEBUG@")
        if(EXISTS "${tbb}")
            get_filename_component(tbb_name "${tbb}" NAME)
            list(APPEND tbb_names "${tbb_name}")
        endif()
    endforeach()
endif()

# Set up a macro to change rpath on given machine
if("@CMAKE_SYSTEM_NAME@" STREQUAL "Darwin")
    macro(change_object_rpath object)
        foreach(directory ${ARGN})
            execute_process(COMMAND
                install_name_tool -add_rpath ${directory} ${object}
                ERROR_QUIET
            )
        endforeach()
    endmacro()
elseif("@CMAKE_SYSTEM_NAME@" STREQUAL "Linux")
    macro(change_object_rpath object)
        set(rpath_string "blob")
        foreach(directory ${ARGN})
            set(rpath_string "${rpath_string}:${directory}")
        endforeach()
        string(REGEX REPLACE "blob:" "" rpath_string "${rpath_string}")
        execute_process(COMMAND
            @PATCHELF_EXECUTABLE@ --set-rpath ${rpath_string} ${object}
        )
    endmacro()
else()
    message(FATAL_ERROR "Don't know how to change RPATH on @CMAKE_SYSTEM_NAME@")
endif()

# Loop over each shared object and add rpaths and change tbbs
foreach(library ${shared_objects})
    get_filename_component(current "${library}" PATH)
    set(current_directories ${directories})
    list(REMOVE_ITEM current_directories "${current}")
    change_object_rpath("${library}" ${current_directories})

    foreach(tbb_name ${tbb_names})
        execute_process(COMMAND
            install_name_tool -change ${tbb_name} @rpath/${tbb_name} ${library}
            ERROR_QUIET
        )
    endforeach()
endforeach()
