file(WRITE ${PROJECT_BINARY_DIR}/original/hello "hello")
file(WRITE ${PROJECT_BINARY_DIR}/original/world "world")
file(WRITE ${PROJECT_BINARY_DIR}/original/something.c "world")
file(WRITE ${PROJECT_BINARY_DIR}/original/andthis.c "world")
file(WRITE ${PROJECT_BINARY_DIR}/original/notthis.h "world")
file(WRITE ${PROJECT_BINARY_DIR}/original/dirA/notthis.c "world")
file(WRITE ${PROJECT_BINARY_DIR}/original/dirA/butthis.h "world")
file(WRITE ${PROJECT_BINARY_DIR}/original/dirA/dirC/andthis.h "world")
file(WRITE ${PROJECT_BINARY_DIR}/original/dirA/dirC/orthis.hh "world")

include(TargetCopyFiles)
include(TestCMake)

add_custom_target(copy ALL)
add_copy_files(copy
    DESTINATION ${PROJECT_BINARY_DIR}/here
    ${PROJECT_BINARY_DIR}/original/hello
    ${PROJECT_BINARY_DIR}/original/world
    VERBOSE
)
add_copy_files(copy
    GLOB ${PROJECT_BINARY_DIR}/original/*.c
    DESTINATION ${PROJECT_BINARY_DIR}/here
    VERBOSE
)
add_copy_directory(copy "${PROJECT_BINARY_DIR}/original"
    RELATIVE "${PROJECT_BINARY_DIR}"
    DESTINATION "${PROJECT_BINARY_DIR}/dir_copies"
    GLOB *.h *.c
    EXCLUDE "notthis.*"
    VERBOSE
)

cmake_test(copy_files)
cmake_test(copy_directories)
