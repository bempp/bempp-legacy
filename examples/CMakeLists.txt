include_directories(${CMAKE_BINARY_DIR}/include)
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/include)
include_directories ("${CMAKE_SOURCE_DIR}/lib")

# Examples
add_executable(tutorial_dirichlet tutorial_dirichlet.cpp)
add_executable(grid_refinement grid_refinement.cpp)
target_link_libraries(tutorial_dirichlet bempp)
target_link_libraries(grid_refinement bempp)
set_target_properties(tutorial_dirichlet PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/bempp/lib")
install(TARGETS tutorial_dirichlet RUNTIME DESTINATION bempp/examples)

set_target_properties(grid_refinement PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/bempp/lib")
install(TARGETS grid_refinement RUNTIME DESTINATION bempp/examples)


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

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/meshes
        DESTINATION bempp/examples)
