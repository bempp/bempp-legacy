# Figures out where to install python stuff
add_subdirectory(doc)
add_subdirectory(bempp)

if (WITH_TESTS)
    add_subdirectory(tests)
endif (WITH_TESTS)
