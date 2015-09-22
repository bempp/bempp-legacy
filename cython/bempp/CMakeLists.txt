include(PythonModule)
include(FilterList)

include_directories("${PYTHON_BINARY_DIR}/bempp/include")

# An umbrella target for cython headeres
add_custom_target(cython-headers)

# Add empty __init__.pxd file so cython considers the include directory
file(WRITE "${PYTHON_BINARY_DIR}/bempp/include/bempp/__init__.pxd")

# Create an __init__.py file that will not be installed later on.
file(WRITE "${PYTHON_BINARY_DIR}/bempp/__init__.py"
    "from pkgutil import extend_path\n__path__ = extend_path(__path__, __name__)")

install_python(FILES __init__.pxd DESTINATION bempp/include/bempp)


#add_python_module(bempp __init__.py 
#    CPP LIBRARIES libbempp)
#add_dependencies(bempp cython-headers)

# Sets path variables for in-build system
set(BEMPP_MESHES "${PROJECT_SOURCE_DIR}/meshes")
set(BEMPP_CMAKE_PATH "${PROJECT_BINARY_DIR}")
configure_file(config.in.py "${PYTHON_BINARY_DIR}/bempp/config.py" @ONLY)

# Set paths for install-tree
function(make_absolute OUTPUT input)
    if("${input}" MATCHES "^/")
        set(output ${input})
    elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
        set(output "/usr/local/${input}")
    elseif("${CMAKE_INSTALL_PREFIX}" MATCHES "^/" )
        set(output "${CMAKE_INSTALL_PREFIX}/${input}")
    else()
        set(output "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/${input}")
    endif()
    set(${OUTPUT} "${output}" PARENT_SCOPE)
endfunction()
make_absolute(BEMPP_MESHES "${SHARE_INSTALL_PATH}/bempp/meshes")
make_absolute(BEMPP_CMAKE_PATH "${SHARE_INSTALL_PATH}/cmake/Bempp")
install_python(FILES "${PYTHON_BINARY_DIR}/bempp/config.py" DESTINATION bempp)

add_subdirectory(core)

