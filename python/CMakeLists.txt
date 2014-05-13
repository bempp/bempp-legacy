# Figures out where to install python stuff
include(PythonInstall)
if(SWIG_FOUND) # Not going anywhere until swig is installed
    add_subdirectory(doc)
    add_subdirectory(bempp)
    add_subdirectory(examples)

    if (WITH_TESTS)
      add_subdirectory(tests)
    endif (WITH_TESTS)
endif()
