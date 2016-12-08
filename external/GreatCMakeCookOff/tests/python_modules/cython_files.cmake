# Check location of binaries in build
set(origin "${PROJECT_BINARY_DIR}/../cython_build")
get_filename_component(origin "${origin}" ABSOLUTE)

set(paths_exist
    "${origin}/python_binary/extension/structure.so"
    "${origin}/python_binary/extension/__init__.py"
)
foreach(pathname ${paths_exist})
    if(NOT EXISTS "${pathname}")
        message(FATAL_ERROR "Path ${pathname} not in build")
    endif()
endforeach()

# Check location of installs
set(paths_exist
    "${origin}/python_install/extension/structure.so"
    "${origin}/python_install/extension/__init__.py"
    "${origin}/python_install/extension/include/extension/structure.pxd"
)
foreach(pathname ${paths_exist})
    if(NOT EXISTS "${pathname}")
        message(FATAL_ERROR "Path ${pathname} not in install")
    endif()
endforeach()

