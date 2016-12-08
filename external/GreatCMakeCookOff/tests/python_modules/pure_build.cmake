find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()
include(PythonModule)

# install and build paths for fake projects
set(PYTHON_BINARY_DIR "${PROJECT_BINARY_DIR}/python_binary"
    CACHE PATH "" FORCE)
set(PYTHON_PKG_DIR "${PROJECT_BINARY_DIR}/python_install"
    CACHE PATH "" FORCE)

# Create fake sources first
if(NOT EXISTS "${PROJECT_BINARY_DIR}/pure")
    file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/pure")
    file(WRITE "${PROJECT_BINARY_DIR}/pure/__init__.py"
        "# Fake dummy package\n"
        "import that\n"
    )
    file(WRITE "${PROJECT_BINARY_DIR}/pure/that.py"
        "# Fake dummy package\n"
        "hello = 'world'\n"
    )

    file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/tests")
    file(WRITE "${PROJECT_BINARY_DIR}/pure/tests/this.py"
        "# Fake dummy package\n"
        "meaning = 42\n"
    )

    file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/noinstall")
    file(WRITE "${PROJECT_BINARY_DIR}/pure/noinstall/__init__.py"
        "# Fake dummy package\n"
        "life = 42.2\n"
    )
endif()

add_python_module("pure"
	${PROJECT_BINARY_DIR}/pure/__init__.py
	${PROJECT_BINARY_DIR}/pure/that.py
)
add_python_module("pure.tests"
    ${PROJECT_BINARY_DIR}/pure/tests/this.py FAKE_INIT)
add_python_module("pure.noinstall"
    "${PROJECT_BINARY_DIR}/pure/noinstall/__init__.py"
    TARGETNAME nope
    NOINSTALL
)

foreach(target pure pure.tests nope)
    if(NOT TARGET ${target})
        message(FATAL_ERROR "Target ${target} does not exist")
    endif()
endforeach()
