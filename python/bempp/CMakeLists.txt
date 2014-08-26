include(PythonModule)
include(MakoFiles)

get_filename_component(header_directory "${CMAKE_CURRENT_SOURCE_DIR}" PATH)
include_directories("${header_directory}")

add_to_python_path("${PROJECT_SOURCE_DIR}/python/mako")
mako_files(options.mako.pyx options.mako.pxd
    OUTPUT_FILES makoed
    DEPENDS "${PROJECT_SOURCE_DIR}/python/mako/aca_options.py"
)
add_python_module(bempp __init__.py ${makoed} EXCLUDE config.in.py CPP)
# CMake does not always pick up dependencies from one add_custom_command to
# another. So we add the dependency explicitly
add_custom_target(bempp-mako DEPENDS ${makoed})
add_dependencies(bempp bempp-mako)

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


add_subdirectory(grid)
add_subdirectory(fiber)
add_subdirectory(assembly)
