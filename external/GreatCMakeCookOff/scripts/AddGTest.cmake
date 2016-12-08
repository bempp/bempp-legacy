# CMake arguments for gtest.
if(NOT MINGW)
    find_package(Threads)
endif()

# sets EXTERNAL_ROOT if not already done
include(PackageLookup)
lookup_package(GTest NOFIND KEEP)

if(PROJECT_USES_CPP11 AND NOT MSVC)
    add_definitions(-DGTEST_LANG_CXX11)
endif(PROJECT_USES_CPP11 AND NOT MSVC)

macro(add_gtest name source)
    include_directories(${EXTERNAL_ROOT}/include)
    # Better, but only works on CMake 2.8.6?
    # get_target_property(THISTEST_INCLUDE test_${name} INCLUDE_DIRECTORIES)
    # set_target_properties(test_${name} PROPERTIES INCLUDE_DIRECTORIES
    #                       "${source_dir}/include;${THISTEST_INCLUDE}")

    add_executable(test_${name} ${source})
    if(EXISTS ${EXTERNAL_ROOT}/lib/gtest.lib)
      target_link_libraries(test_${name} ${EXTERNAL_ROOT}/lib/gtest.lib)
    elseif(EXISTS ${EXTERNAL_ROOT}/lib/libgtest.a)
      target_link_libraries(test_${name} ${EXTERNAL_ROOT}/lib/libgtest.a)
    elseif(NOT GTest_LOOKUP_BUILD)
      message(FATAL_ERROR "Could not find gtest")
    endif()
    if(CMAKE_THREAD_LIBS_INIT)
      target_link_libraries(test_${name} ${CMAKE_THREAD_LIBS_INIT})
    endif(CMAKE_THREAD_LIBS_INIT)

    add_dependencies(test_${name} Lookup-GTest)
    if(NOT "${ARGN}" STREQUAL "")
      target_link_libraries(test_${name} ${ARGN})
    endif(NOT "${ARGN}" STREQUAL "")

    add_test(cxx_${name} ${EXECUTABLE_OUTPUT_PATH}/test_${name}
                --gtest_output=xml:${CMAKE_BINARY_DIR}/test-results/test_${name}.xml)
    set_tests_properties(cxx_${name} PROPERTIES LABELS "gtest")
endmacro()
