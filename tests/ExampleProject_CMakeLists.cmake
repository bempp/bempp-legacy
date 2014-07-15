# Test template Bempp project
cmake_minimum_required(VERSION 2.8)
project(tut_test CXX)

find_package(Bempp @Bempp_VERSION@ EXACT REQUIRED)

# Tests that find_package(Bempp) declares what we think it should
# Delete the next few lines if using this as BEM++ project template
if(PROJECT_INCLUSION_UNIT_TESTS)
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
    include(ProjectInclusionUnitTests)
endif()

# Add an executable that links to BEM++
include_directories("${BEMPP_INCLUDE_DIRS}")
add_executable(tut_test
    "@PROJECT_SOURCE_DIR@/examples/cpp/tutorial_dirichlet.cpp")
target_link_libraries(tut_test bempp)
