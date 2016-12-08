find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()
include(AddPyTest)
enable_testing()

# install and build paths for fake projects
set(PYTHON_BINARY_DIR "${PROJECT_BINARY_DIR}/python_binary"
    CACHE PATH "" FORCE)
set(PYTHON_PKG_DIR "${PROJECT_BINARY_DIR}/python_install"
    CACHE PATH "" FORCE)
set(LOCAL_PYTHON_EXECUTABLE "@PROJECT_BINARY_DIR@/localpython.sh")

# Create fake sources first
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/hackage")
    file(MAKE_DIRECTORY "${PROJECT_SOURCE_DIR}/hackage")
endif()
file(WRITE "${PROJECT_SOURCE_DIR}/hackage/test_this.py"
    "def test_something():\n"
    "   assert True\n"
    "def test_that():\n"
    "   assert 1 == 1"
)
file(WRITE "${PROJECT_SOURCE_DIR}/hackage/test_that.py"
    "def test_that():\n"
    "   assert False"
)

# Command-line arguments
file(WRITE "${PROJECT_SOURCE_DIR}/hackage/conftest.py"
    "from py.test import fixture\n"
    "def pytest_addoption(parser):\n"
    "    parser.addoption('--cmdl', action='store')\n\n\n"
    "@fixture\n"
    "def cmdl(request):\n"
    "    return request.config.getoption('--cmdl')\n"
)
file(WRITE "${PROJECT_SOURCE_DIR}/hackage/test_cmdl.py"
    "def test_that(cmdl):\n"
    "   assert cmdl == 'an option'"
)

# cython tests
file(WRITE "${PROJECT_SOURCE_DIR}/hackage/_cython.pyx"
    "def test_something():\n"
    "   cdef int value = 5\n"
    "   assert value == 5\n"
    "   assert value != 6\n"
)
file(WRITE "${PROJECT_SOURCE_DIR}/hackage/_cython_cpp.pyx"
    "from libcpp.vector cimport vector\n"
    "def test_cpp():\n"
    "   cdef vector[int] *value = new vector[int](1, 3)\n"
    "   assert value.at(0) == 3\n"
    "   assert value.size() == 1\n"
    "   del value"
)
file(WRITE "${PROJECT_SOURCE_DIR}/hackage/test_cython.py"
    "def test_something():\n"
    "   from _cython import test_something\n"
    "   test_something()\n"
    "def test_cpp():\n"
    "   from _cython_cpp import test_cpp\n"
    "   test_cpp()\n"
)

find_package(CoherentPython REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

setup_pytest("@EXTERNAL_ROOT@/python" "@PROJECT_BINARY_DIR@/py.test.sh")

add_pytest(hackage/test_*.py hackage/*.pyx
    EXCLUDE hackage/test_cmdl.py
    PREFIX "hackage"
    CPP
)
add_pytest(hackage/test_cmdl.py hackage/conftest.py
    PREFIX "hackage" CMDLINE "--cmdl=an option")
add_pytest(hackage/test_cmdl.py hackage/conftest.py
    PREFIX "hackage.fails" CMDLINE "--cmdl=nope")
