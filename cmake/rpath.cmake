# Sets RPATH so that the executable
# and python extensions can find libpurify.dylib
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    if(NOT DEFINED CMAKE_MACOSX_RPATH)
        set(CMAKE_MACOSX_RPATH TRUE CACHE BOOL
            "Use rpath to make libraries relocatable")
    endif()
endif()

# Set RPATH to location where libraries will be installed,
# unless it is already part of the platform
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
    "${LIBRARY_INSTALL_PATH}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    list(APPEND CMAKE_INSTALL_RPATH "${LIBRARY_INSTALL_PATH}")
    list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)
    set(CMAKE_INSTALL_RPATH
        ${CMAKE_INSTALL_RPATH}
        CACHE PATH
        "Path where shared libraries should be installed"
        FORCE
    )
endif()
