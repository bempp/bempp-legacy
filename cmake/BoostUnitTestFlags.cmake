# Checks that boost whether to have BOOST_TEST_DYN_LINK defined or not.
# It should be defined if linking to a dynamic boost unit test library,
# and absent when linking to a static boost unit test library.
if(NOT Boost_FOUND OR DEFINED Boost_UNIT_TEST_CXXFLAGS)
    return()
endif()
set(rootdir "${PROJECT_BINARY_DIR}/CMakeFiles/boostStuff")
file(WRITE "${rootdir}/main.cc"
    "#define BOOST_TEST_MODULE something\n"
    "#include <boost/test/unit_test.hpp>\n"
)
file(WRITE "${rootdir}/CMakeLists.txt"
    "cmake_minimum_required(VERSION 2.8)\n"
    "project(boostflags CXX)\n"
    "set(Boost_INCLUDE_DIR \"${Boost_INCLUDE_DIR}\")\n"
    "set(Boost_LIBRARIES \"${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}\")\n"
    "include_directories(\${Boost_INCLUDE_DIR})\n"
    "if(DO_DYNAMIC)\n"
    "   add_definitions(-DBOOST_TEST_DYN_LINK)\n"
    "endif()\n"
    "add_executable(flagme main.cc)\n"
    "target_link_libraries(flagme \${Boost_LIBRARIES})\n"
)
try_compile(without_flag
    "${rootdir}/without" "${rootdir}"
    boostflags flagme
    CMAKE_FLAGS -DDO_DYNAMIC=FALSE
    OUTPUT_VARIABLE OUTPUT_WITHOUT
)
try_compile(with_flag
    "${rootdir}/with" "${rootdir}"
    boostflags flagme
    CMAKE_FLAGS -DDO_DYNAMIC=TRUE
    OUTPUT_VARIABLE OUTPUT_WITH
)
if(without_flag AND NOT with_flag)
    set(result FALSE)
elseif(NOT without_flag AND with_flag)
    set(result -DBOOST_TEST_DYN_LINK)
else()
    message("without: ${OUTPUT_WITHOUT}")
    message("with: ${OUTPUT_WITH}")
    message("${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")
    message(FATAL_ERROR "Could not determine boost unit test flags")
endif()
set(Boost_UNIT_TEST_CXXFLAGS ${result}
    CACHE STRING "Compile flags for boost unit tests")
