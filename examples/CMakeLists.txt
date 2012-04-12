include_directories ("${PROJECT_SOURCE_DIR}/lib")
include_directories(${CMAKE_BINARY_DIR}/include/common)

if (WITH_TRILINOS)
    include_directories(${TRILINOS_INCLUDE_DIR})
    add_executable(trilinos trilinos.cpp meshes.cpp)
    target_link_libraries(trilinos bempp
        ${BLAS_LIBRARIES} # must be put before the Trilinos libraries,
                          # otherwise BLAS behaves incorrectly
        ${STRATIMIKOS_LIBRARY}
        ${STRATIMIKOS_BELOS_LIBRARY}
        ${HDF5_LIBRARY})
endif ()

if (WITH_ALUGRID)
    add_definitions(-DWITH_ALUGRID -DHAVE_ALUGRID=1 -DENABLE_ALUGRID)
    include_directories(${ALUGRID_INCLUDE_DIRS})
    add_executable(surface_mesh_extraction surface_mesh_extraction.cpp)
    target_link_libraries(surface_mesh_extraction bempp ${ALUGRID_LIBRARY})
endif ()

add_executable(traversal traversal.cpp)
add_executable(traversal_benchmark traversal_benchmark.cpp)
add_executable(gmsh_import gmsh_import.cpp)
add_executable(vtk_export vtk_export.cpp)
add_executable(assembly assembly.cpp meshes.cpp)
add_executable(rhs_assembly rhs_assembly.cpp)
add_executable(superposition superposition.cpp meshes.cpp)
add_executable(calderon calderon.cpp meshes.cpp)
add_executable(physical_test physical_test.cpp meshes.cpp)
add_executable(tutorial_dirichlet tutorial_dirichlet.cpp)
target_link_libraries(traversal bempp)
target_link_libraries(traversal_benchmark bempp)
target_link_libraries(gmsh_import bempp)
target_link_libraries(vtk_export bempp)
target_link_libraries(assembly bempp)
target_link_libraries(rhs_assembly bempp)
target_link_libraries(superposition bempp)
target_link_libraries(calderon bempp)
target_link_libraries(physical_test bempp)
target_link_libraries(tutorial_dirichlet bempp)

if (WITH_ALUGRID)
    add_executable(extract_surface_grid extract_surface_grid.cpp)
    target_link_libraries (extract_surface_grid bempp)
endif ()

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
