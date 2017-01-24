foreach(filename something.c andthis.c dirA/butthis.h dirA/dirC/andthis.h)
    if(NOT EXISTS @PROJECT_BINARY_DIR@/dir_copies/original/${filename})
        message(FATAL_ERROR "Did not copy original/original/${filename}")
    endif()
endforeach()

foreach(filename hello world notthis.h  dirA/notthis.c dirA/dirC/orthis.hh)
    if(EXISTS @PROJECT_BINARY_DIR@/dir_copies/original/${filename})
        message(FATAL_ERROR "Should not have copied original/${filename}")
    endif()
endforeach()

