include_directories(${CMAKE_BINARY_DIR}/include)
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/include)
include_directories ("${CMAKE_SOURCE_DIR}/lib")

add_executable(tutorial_dirichlet tutorial_dirichlet.cpp)
target_link_libraries(tutorial_dirichlet bempp)

if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    file(GLOB_RECURSE EXAMPLE_MESHES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
        *.msh *.gmsh sphere.txt)
    foreach (f ${EXAMPLE_MESHES})
        add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${f}
            COMMAND ${CMAKE_COMMAND}
                copy -E
                    ${CMAKE_CURRENT_SOURCE_DIR}/${f}
                    ${CMAKE_CURRENT_BINARY_DIR}/${f}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${f}
            # TODO: remove the project's root directory from the printed path.
            COMMENT "Copying file ${CMAKE_CURRENT_SOURCE_DIR}/${f}"
            VERBATIM
        )
    endforeach ()

    foreach (f ${EXAMPLE_MESHES})
        set(EXAMPLE_MESH_DEPS
            "${EXAMPLE_MESH_DEPS};${CMAKE_CURRENT_BINARY_DIR}/${f}")
    endforeach ()
    add_custom_target(copy_example_deps ALL
        DEPENDS ${EXAMPLE_MESH_DEPS})
endif ()
