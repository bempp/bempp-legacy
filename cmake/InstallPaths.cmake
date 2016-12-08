# Sets up install paths
# Mostly, add option to easily install as a python package
if("${PYTHON_PKG_DIR}" STREQUAL "")
    # Should not happen if PythonInstall has already been included
    message(FATAL_ERROR "PYTHON_PKG_DIR is not set")
endif()
set(prefix "${PYTHON_PKG_DIR}/bempp/")
set(LIBRARY_INSTALL_PATH "${prefix}lib")
set(INCLUDE_INSTALL_PATH "${prefix}include")
set(SHARE_INSTALL_PATH "${prefix}share")
set(RUNTIME_INSTALL_PATH "${prefix}bin")
set(DOC_INSTALL_PATH "${prefix}doc")
