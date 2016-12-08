cmake_test(checkisnan SOURCE)
cmake_test(checkcpp11flags)
cmake_test(local_include)
cmake_test(self_lookup)
cmake_test(cached_variables)
cmake_test(delete_variables)
cmake_test(passon_variables)
cmake_test(python_package)
cmake_test(add_to_path)
cmake_test(create_environment_script)
cmake_test(output_filename)
cmake_test(configure_files)
if(CMAKE_VERSION VERSION_GREATER "3.0" OR CMAKE_VERSION VERSION_EQUAL "3.0")
  cmake_test(hunter HEADER_LINES
    "find_package(GreatCMakeCookOff NO_MODULE NO_DEFAULT_PATH PATHS \${cookoff_path} REQUIRED)\n"
    "initialize_cookoff()\n"
    "include(AddHunter)\n"
    "initialize_hunter(HUNTER_ROOT \${CMAKE_BINARY_DIR}/hunter_root)\n"
  )
endif()
if(TEST_FIND_MKL)
  cmake_test(find_mkl SOURCE NOEXEC)
  cmake_test(find_mkl_scalapack SOURCE NOEXEC)
endif()

cmake_test(lookup_eigen NOEXEC)
find_program(Hg_EXECUTABLE hg)
if("${CMAKE_VERSION}" VERSION_GREATER 2.8.9 AND Hg_EXECUTABLE)
    cmake_test(use_external)
    cmake_test(use_external_download_by_default)
    set_tests_properties(
        cmake_test_use_external cmake_test_use_external_download_by_default
        PROPERTIES DEPENDS cmake_test_lookup_eigen
    )
endif()
find_program(svn_EXECUTABLE svn)
if(svn_EXECUTABLE)
    cmake_test(addgtest NOEXEC SOURCE
        --test-command ${CMAKE_MAKE_PROGRAM} test)
endif()
cmake_test(addcatchtest NOEXEC SOURCE --test-command ${CMAKE_MAKE_PROGRAM} test)


add_subdirectory(cpp11)
add_subdirectory(copy_files)
add_subdirectory(python_modules)
add_subdirectory(py.test)
add_subdirectory(mako)
