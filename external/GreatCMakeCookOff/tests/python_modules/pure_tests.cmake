# Check location of binaries in build
set(origin "${PROJECT_BINARY_DIR}/../pure_build")
set(paths_exist
    "${origin}/python_binary/pure"
    "${origin}/python_binary/pure/__init__.py"
    "${origin}/python_binary/pure/that.py"
    "${origin}/python_install/pure/tests/"
    "${origin}/python_binary/pure/tests/__init__.py"
    "${origin}/python_binary/pure/tests/this.py"
    "${origin}/python_binary/pure/noinstall/__init__.py"
)
foreach(pathname ${paths_exist})
    if(NOT EXISTS "${pathname}")
        message(FATAL_ERROR "Path ${pathname} not in build")
    endif()
endforeach()

# Check location of installs
set(paths_exist
    "${origin}/python_install/pure"
    "${origin}/python_install/pure/__init__.py"
    "${origin}/python_install/pure/that.py"
    "${origin}/python_install/pure/tests/"
    "${origin}/python_install/pure/tests/__init__.py"
    "${origin}/python_install/pure/tests/this.py"
)
foreach(pathname ${paths_exist})
    if(NOT EXISTS "${pathname}")
        message(FATAL_ERROR "Path ${pathname} not in install")
    endif()
endforeach()

# Check location does not exist
set(paths_exist
    "${origin}/python_install/noinstall"
    "${origin}/python_install/noinstall/__init__.py"
)
foreach(pathname ${paths_exist})
    if(EXISTS "${pathname}")
        message(FATAL_ERROR "Path ${pathname} should not exist")
    endif()
endforeach()
