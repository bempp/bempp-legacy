include_directories ("${PROJECT_SOURCE_DIR}/lib")
add_executable(example1 example1.cpp)
target_link_libraries(example1 bempp)
