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
    endif()
endforeach()
list(REMOVE_DUPLICATES shared_objects)
list(REMOVE_DUPLICATES directories)

# Gets list of tbbs
unset(tbb_names)
foreach(tbb "@TBB_LIBRARY@" "@TBB_LIBRARY_DEBUG@"
    "@TBB_MALLOC_LIBRARY@" "@TBB_MALLOC_LIBRARY_DEBUG@")
    if(EXISTS "${tbb}")
        get_filename_component(tbb_name "${tbb}" NAME)
        list(APPEND tbb_names "${tbb}")
    endif()
endforeach()
# Loop over each shared object and add rpaths and change tbbs
foreach(library ${shared_objects})
    get_filename_component(current "${library}" PATH)
    set(current_directories ${directories})
    list(REMOVE_ITEM current_directories "${current}")
    foreach(directory ${current_directories})
        file(RELATIVE_PATH relative "${current}" "${directory}/fake")
        get_filename_component(relative "${relative}" PATH)
        execute_process(COMMAND
            install_name_tool -add_rpath ${relative} ${library}
            ERROR_QUIET
        )
    endforeach()

    foreach(tbb_name ${tbb_names})
        execute_process(COMMAND
            install_name_tool -change ${tbb_name} @rpath/${tbb_name} ${library}
            ERROR_QUIET
        )
    endforeach()
endforeach()
