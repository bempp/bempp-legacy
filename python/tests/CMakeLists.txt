add_subdirectory(unit)
if (WITH_INTEGRATION_TESTS)
   add_subdirectory(integration)
endif ()
