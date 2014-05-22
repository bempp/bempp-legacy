set(tbbs "@TBB_LIBRARY@" "@TBB_LIBRARY_DEBUG@" "@TBB_MALLOC_LIBRARY@"
    "@TBB_MALLOC_LIBRARY_DEBUG@")
file(GLOB_RECURSE libraries "*.so" "*.dylib")
foreach(tbb ${tbbs})
    if(NOT "${tbb}" STREQUAL "" AND EXISTS "${tbb}")
        get_filename_component(tbb_name "${tbb}" NAME)
        foreach(library ${libraries})
            execute_process(COMMAND
                install_name_tool -change ${tbb_name} @rpath/${tbb_name} ${library}
            )
        endforeach()
    endif()
endforeach()
