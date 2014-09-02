include(PythonModule)
include(MakoFiles)
include(FilterList)

include_directories("${PYTHON_BINARY_DIR}/bempp/include")

add_to_python_path("${PROJECT_SOURCE_DIR}/python/mako")
set(makoes options.mako.pyx options.mako.pxd)
mako_files(${makoes}
    OUTPUT_FILES makoed
    DEPENDS "${PROJECT_SOURCE_DIR}/python/mako/bempp_options.py"
    TARGETNAME bempp-mako
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp"
)

# An umbrella target for cython headeres
add_custom_target(cython-headers)
add_dependencies(cython-headers bempp-mako)

# Add empty __init__.pxd file so cython considers the include directory
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/__init__.pxd")
list(APPEND makoed "${CMAKE_CURRENT_BINARY_DIR}/__init__.pxd")

# headers are setup in python package first so that they can be accessed by
# cython. The main issue is for mako-generated headers.
split_list(pxd_headers pyx_sources makoed ".*\\.pxd")
install_python(FILES ${headers} DESTINATION bempp/include/bempp)


add_python_module(bempp __init__.py ${pyx_sources} CPP LIBRARIES libbempp)
add_dependencies(bempp cython-headers)

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
configure_file(config.in.py "${CMAKE_CURRENT_BINARY_DIR}/config.py" @ONLY)
install_python(FILES "${CMAKE_CURRENT_BINARY_DIR}/config.py" DESTINATION bempp)


add_subdirectory(utils)
add_subdirectory(grid)
add_subdirectory(fiber)
add_subdirectory(space)
