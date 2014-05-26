# Examples
add_executable(tutorial_dirichlet tutorial_dirichlet.cpp)
target_link_libraries(tutorial_dirichlet bempp)

# Meshes
file(GLOB_RECURSE EXAMPLE_MESHES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    *.msh sphere.txt)
if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    foreach (f ${EXAMPLE_MESHES})
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${f}
                       ${CMAKE_CURRENT_BINARY_DIR}/${f}
                       COPYONLY)
    endforeach ()
endif ()

install(TARGETS tutorial_dirichlet
    EXPORT BemppTargets
    RUNTIME
    DESTINATION ${RUNTIME_INSTALL_PATH}/bempp/examples)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/meshes
    DESTINATION ${SHARE_INSTALL_PATH}/bempp/examples)
