cmake_test(mako_build NOEXEC --build-target install)
cmake_test(mako_tests)
set_tests_properties(cmake_test_mako_build cmake_test_mako_tests
    PROPERTIES LABELS "mako")
set_tests_properties(cmake_test_mako_tests PROPERTIES
    DEPENDS cmake_test_mako_build)
