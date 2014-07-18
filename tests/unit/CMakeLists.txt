include(BoostUnitTestFlags)
if(Boost_UNIT_TEST_CXXFLAGS)
    add_definitions(${Boost_UNIT_TEST_CXXFLAGS})
endif()

# Adds a bash script to run cpp unit tests
create_environment_script(
    WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/tests/unit"
    PATH "${PROJECT_BINARY_DIR}/run_cpptests.sh"
)


# builds the command to link to boost::unittest
set(btest_libs)
if(Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG)
    list(APPEND btest_libs debug ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG})
endif()
if(Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE)
    list(APPEND btest_libs optimized
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})
elseif(Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG)
    list(REMOVE_AT btest_libs 0)
endif()
if("${btest_libs}" STREQUAL "")
    message(FATAL_ERROR "No boost unit test library")
endif()

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

# tests that need extra sources should be declared below
function(extra_sources SOURCE)
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
        list(APPEND extras "assembly/create_regular_grid.cpp")
    endif()
    if("${filename}" STREQUAL "default_direct_solver"
            OR "${filename}" STREQUAL "default_iterative_solver")
        list(APPEND extras "linalg/laplace_3d_dirichlet_fixture.cpp")
    endif()
    if("${filename}" STREQUAL "entity"
        OR "${filename}" STREQUAL "grid"
        OR "${filename}" STREQUAL "grid_view"
        OR "${filename}" STREQUAL "geometry"
        OR "${filename}" STREQUAL "id_set"
        OR "${filename}" STREQUAL "grid_factory"
        OR "${filename}" STREQUAL "index_set"
    )
        list(APPEND extras "grid/simple_triangular_grid_manager.cpp")
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
    extra_sources("${source}")
    add_executable(${testname} main.cpp "${source}" ${extras})
    depends_on_lookups(${testname})
    target_link_libraries(${testname} bempp)
    target_link_libraries(${testname} ${btest_libs})
    
    add_test(NAME ${testname} COMMAND ${test_sh} $<TARGET_FILE:${testname}>)
    set_tests_properties(${testname} PROPERTIES LABELS "cpp;unit;${basename}")
endforeach()

# Meshes
file(GLOB_RECURSE TEST_MESHES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    *.msh)
if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    foreach (f ${TEST_MESHES})
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${f}
                       ${CMAKE_CURRENT_BINARY_DIR}/${f}
                       COPYONLY)
    endforeach ()
endif ()
