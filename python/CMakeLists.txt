add_subdirectory(doc)
add_subdirectory(bempp)
add_subdirectory(examples)

if (WITH_TESTS)
  add_subdirectory(tests)
endif (WITH_TESTS)
