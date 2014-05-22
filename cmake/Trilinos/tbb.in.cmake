set(tbbs "@TBB_LIBRARY@" "@TBB_LIBRARY_DEBUG@" "@TBB_MALLOC_LIBRARY@"
    "@TBB_MALLOC_LIBRARY_DEBUG@")
file(GLOB_RECURSE dylibs "*.dylib")
file(GLOB_RECURSE pyexts "*.so")
unset(actual_tbbs)
foreach(tbb ${tbbs})
    if(NOT "${tbb}" STREQUAL "" AND EXISTS "${tbb}")
        list(APPEND actual_tbbs "${tbb}")
    endif()
endforeach()
foreach(library ${dylibs} ${pyexts})
    foreach(tbb ${actual_tbbs})
        get_filename_component(tbb_name "${tbb}" NAME)
        execute_process(COMMAND
            install_name_tool -change ${tbb_name} @rpath/${tbb_name} ${library}
        )
    endforeach()
endforeach()
