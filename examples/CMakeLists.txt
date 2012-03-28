include_directories ("${PROJECT_SOURCE_DIR}/lib")
add_executable(traversal traversal.cpp)
add_executable(traversal_benchmark traversal_benchmark.cpp)
add_executable(gmsh_import gmsh_import.cpp)
add_executable(vtk_export vtk_export.cpp)
add_executable(assembly assembly.cpp meshes.cpp)
add_executable(rhs_assembly rhs_assembly.cpp)
#add_executable(show_shape_functions show_shape_functions.cpp)
target_link_libraries(traversal bempp)
target_link_libraries(traversal_benchmark bempp)
target_link_libraries(gmsh_import bempp)
target_link_libraries(vtk_export bempp)
target_link_libraries(assembly bempp)
target_link_libraries(rhs_assembly bempp)
#target_link_libraries(show_shape_functions bempp)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/head.gmsh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/simple_mesh_2_elements.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/simple_mesh_9_elements.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/two_disjoint_triangles.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    file(GLOB_RECURSE EXAMPLE_MESHES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.msh)
    foreach (f ${EXAMPLE_MESHES})
        add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${f}
            COMMAND ${CMAKE_COMMAND} copy -E ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${CMAKE_CURRENT_BINARY_DIR}/${f}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${f}
            # TODO: remove the project's root directory from the printed path.
            COMMENT "Copying file ${CMAKE_CURRENT_SOURCE_DIR}/${f}"
            VERBATIM
        )
    endforeach ()

    foreach (f ${EXAMPLE_MESHES})
        set(EXAMPLE_MESH_DEPS "${EXAMPLE_MESH_DEPS};${CMAKE_CURRENT_BINARY_DIR}/${f}")
    endforeach ()
    add_custom_target(copy_example_meshes ALL DEPENDS ${EXAMPLE_MESH_DEPS})
endif ()
