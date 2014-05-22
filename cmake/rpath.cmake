# Sets RPATH so that the executable
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    if(NOT DEFINED CMAKE_MACOSX_RPATH)
        set(CMAKE_MACOSX_RPATH TRUE CACHE BOOL
            "Use rpath to make libraries relocatable")
    endif()
endif()

function(add_to_rpath)
    foreach(path ${ARGN})
        list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${path}" system_dir)
        if("${system_dir}" STREQUAL "-1")
            list(APPEND CMAKE_INSTALL_RPATH "${path}")
        endif()
    endforeach()

    list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)
    set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} PARENT_SCOPE)
endfunction()

function(change_tbb_install_name target)
    if(NOT "${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin" OR NOT TBB_LIBRARY)
        return()
    endif()
    unset(commands)
    foreach(version "" _DEBUG)
        foreach(library "" MALLOC_)
            set(varname TBB_${library}LIBRARY${version})
            if(${varname})
                get_filename_component(library "${${varname}}" NAME)
                list(APPEND commands COMMAND
                    install_name_tool -change
                    ${library} @rpath/${library} $<TARGET_LINKER_FILE:${target}>
                )
            endif()
        endforeach()
    endforeach()

    add_custom_command(TARGET ${target} POST_BUILD
        ${commands}
        COMMENT "${target} to locate TBB libraries with @rpath"
    )
endfunction()
