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
    set(CMAKE_INSTALL_RPATH
        ${CMAKE_INSTALL_RPATH}
        CACHE PATH
        "Path where shared libraries should be installed"
        FORCE
    )
endfunction()
