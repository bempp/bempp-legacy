include(TargetCopyFiles)
include(BoostUnitTestFlags)
if(Boost_UNIT_TEST_CXXFLAGS)
    add_definitions(${Boost_UNIT_TEST_CXXFLAGS})
endif()

# Adds a bash script to run cpp unit tests
create_environment_script(
    WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/tests/unit"
    PATH "${PROJECT_BINARY_DIR}/run_cpptests.sh"
)



# Create a list of tests
file(GLOB_RECURSE TEST_SOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    test_*.cpp)
if(NOT WITH_AHMED)
    list(REMOVE_ITEM TEST_SOURCES
       assembly/test_aca_mode_assembly.cpp
	   assembly/test_discrete_aca_boundary_operator.cpp
	   assembly/test_integral_evaluation_with_blas.cpp
	   assembly/test_synthetic_boundary_operator.cpp
    )
endif()

# declare fixture libraries that will be used by tests
add_library(grid_fixture STATIC "assembly/create_regular_grid.cpp")
add_library(dirichlet_fixture STATIC "linalg/laplace_3d_dirichlet_fixture.cpp")
add_library(manager_fixture STATIC "grid/simple_triangular_grid_manager.cpp")
foreach(fixture grid dirichlet manager)
    depends_on_lookups(${fixture}_fixture)
    add_dependencies(${fixture}_fixture copy_headers)
endforeach()

# tests that need extra dependencies should be declared below
function(extra_dependencies SOURCE)
    unset(extras)
    get_filename_component(filename "${source}" NAME_WE)
    string(REGEX REPLACE "test_?" "" filename "${filename}")
    if("${filename}" STREQUAL "discrete_aca_boundary_operator"
        OR "${filename}" STREQUAL "discrete_boundary_operator_sum"
        OR "${filename}" STREQUAL "discrete_dense_boundary_operator"
        OR "${filename}" STREQUAL "discrete_inverse_sparse_boundary_operator"
        OR "${filename}" STREQUAL "discrete_null_boundary_operator"
        OR "${filename}" STREQUAL "discrete_sparse_boundary_operator"
        OR "${filename}" STREQUAL "sparse_cholesky"
        OR "${filename}" STREQUAL "raviart_thomas_0_vector_space"
    )
        list(APPEND extras grid_fixture)
    endif()
    if("${filename}" STREQUAL "default_direct_solver"
            OR "${filename}" STREQUAL "default_iterative_solver")
        list(APPEND extras dirichlet_fixture)
    endif()
    if("${filename}" STREQUAL "entity"
        OR "${filename}" STREQUAL "grid"
        OR "${filename}" STREQUAL "grid_view"
        OR "${filename}" STREQUAL "geometry"
        OR "${filename}" STREQUAL "id_set"
        OR "${filename}" STREQUAL "grid_factory"
        OR "${filename}" STREQUAL "index_set"
    )
        list(APPEND extras manager_fixture)
    endif()
    set(extras ${extras} PARENT_SCOPE)
endfunction()

# Creates each test
set(test_sh "${PROJECT_BINARY_DIR}/run_cpptests.sh")
foreach(source ${TEST_SOURCES})
    get_filename_component(filename "${source}" NAME_WE)
    string(REGEX REPLACE "test_?" "" filename "${filename}")
    get_filename_component(basename "${source}" PATH)
    get_filename_component(basename "${basename}" NAME)


    set(testname ${basename}_${filename})
    add_executable(${testname} main.cpp "${source}")
    depends_on_lookups(${testname})
    extra_dependencies("${source}")
    target_link_libraries(${testname} libbempp ${extras})
    if(Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG)
        target_link_libraries(${testname}
        debug ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG}
    )
    endif()
    if(Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE)
        target_link_libraries(${testname}
        optimized ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE}
    )
    endif()

    add_test(NAME ${testname} COMMAND ${test_sh} $<TARGET_FILE:${testname}>)
    set_tests_properties(${testname} PROPERTIES LABELS "cpp;unit;${basename}")
endforeach()

# Meshes
add_custom_target(copy_meshes ALL)
add_copy_directory(copy_meshes "${CMAKE_CURRENT_SOURCE_DIR}"
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
    GLOB *.msh
)
