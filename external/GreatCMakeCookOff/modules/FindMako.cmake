# First check for python executable
if(NOT PYTHON_EXECUTABLE)
  find_program(PYTHON_EXECUTABLE python)
  if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "Please check for python before searching for mako")
  endif()
endif()
# If mako was installed by this process previously, check that it is still there
if(MAKO_EXTERNAL_INSTALL)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} "-c \"import mako\""
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/external/python
    ERROR_VARIABLE DUMMY
    OUTPUT_VARIABLE DUMMY
    RESULT_VARIABLE MAKO_STILL_THERE
  )
  if(MAKO_STILL_THERE)
    return()
  endif()
endif()
# If mako is not found ...
if(NOT FOUND_MAKO)
  # Then look for it by importing it
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import mako"
    ERROR_VARIABLE DUMMY
    OUTPUT_VARIABLE DUMMY
    RESULT_VARIABLE FOUND_MAKO
  )
  if(FOUND_MAKO EQUAL 0)
    message(STATUS "[Mako] Found")
    set(FOUND_MAKO TRUE CACHE INTERNAL "Found mako python package")
    set(MAKO_EXTERNAL_INSTALL FALSE CACHE INTERNAL "Manual install of mako")
    return()
  endif()
else()
  # If mako already found then continue
  # We have checked that it was still there, if we installed it
  return()
endif()

message(STATUS "[Mako] Not found")
# Now try and install mako using pip
if(NOT PIP_EXECUTABLE)
  find_program(PIP_EXECUTABLE pip)
endif()
if(NOT PIP_EXECUTABLE)
  message(FATAL_ERROR "cannot install package without pip")
else()
  message(STATUS "[pip] found: ${PIP_EXECUTABLE}")
endif()
# We install it in the build directory
execute_process(
  COMMAND ${PIP_EXECUTABLE} install mako
           --install-option=--install-purelib=${PROJECT_BINARY_DIR}/external/python
           --install-option=--install-scripts=${PROJECT_BINARY_DIR}/external/python
           --install-option=--prefix=${PROJECT_BINARY_DIR}/external/python
  OUTPUT_VARIABLE PIP_OUTPUT
  ERROR_VARIABLE PIP_ERROR
  RESULT_VARIABLE PIP_INSTALLATION_WORKED
)
if(NOT PIP_INSTALLATION_WORKED EQUAL 0)
  message(STATUS "${PIP_OUTPUT}")
  message(STATUS "${PIP_ERROR}")
  message(FATAL_ERROR "Could not install mako. Please see error message above")
else()
  message(STATUS "[Mako] installed in ${PROJECT_BINARY_DIR}/external/python")
endif()

execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c "import mako"
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/external/python
  ERROR_VARIABLE DUMMY
  OUTPUT_VARIABLE DUMMY
  RESULT_VARIABLE FOUND_MAKO
)
if(FOUND_MAKO EQUAL 0)
  message(STATUS "[Mako] now available in build directory")
  set(FOUND_MAKO TRUE CACHE INTERNAL "Found mako python package")
  set(MAKO_EXTERNAL_INSTALL TRUE CACHE INTERNAL "Manual install of mako")
  return()
else()
  message(FATAL_ERROR "Could not install mako")
endif()
