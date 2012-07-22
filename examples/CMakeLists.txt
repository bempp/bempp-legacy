include_directories(${CMAKE_BINARY_DIR}/include)
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/include)
include_directories ("${CMAKE_SOURCE_DIR}/lib")

add_executable(dirichlet dirichlet.cpp meshes.cpp)
add_executable(neumann neumann.cpp meshes.cpp)
add_executable(dot_neumann dot_neumann.cpp meshes.cpp)
add_executable(dot_two_layers dot_two_layers.cpp meshes.cpp)
target_link_libraries(dirichlet bempp)
target_link_libraries(neumann bempp)
target_link_libraries(dot_neumann bempp)
target_link_libraries(dot_two_layers bempp)

if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    file(GLOB_RECURSE EXAMPLE_MESHES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.msh *gmsh)
    file(GLOB_RECURSE TRILINOS_XMLS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.xml)
    foreach (f ${EXAMPLE_MESHES} ${TRILINOS_XMLS})
        add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${f}
            COMMAND ${CMAKE_COMMAND} copy -E ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${CMAKE_CURRENT_BINARY_DIR}/${f}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${f}
            # TODO: remove the project's root directory from the printed path.
            COMMENT "Copying file ${CMAKE_CURRENT_SOURCE_DIR}/${f}"
            VERBATIM
        )
    endforeach ()

    foreach (f ${EXAMPLE_MESHES} ${TRILINOS_XMLS})
        set(EXAMPLE_MESH_DEPS "${EXAMPLE_MESH_DEPS};${CMAKE_CURRENT_BINARY_DIR}/${f}")
    endforeach ()
    add_custom_target(copy_example_deps ALL DEPENDS ${EXAMPLE_MESH_DEPS} ${TRILINOS_XMLS})
endif ()
