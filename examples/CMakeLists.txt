include_directories ("${PROJECT_SOURCE_DIR}/lib")
add_executable(traversal traversal.cpp)
target_link_libraries(traversal bempp dunecommon dunegrid)

