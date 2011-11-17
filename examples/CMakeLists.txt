include_directories ("${PROJECT_SOURCE_DIR}/lib")
add_executable(traversal traversal.cpp)
add_executable(traversal_benchmark traversal_benchmark.cpp)
add_executable(pure_traversal_benchmark pure_traversal_benchmark.cpp)
target_link_libraries(traversal bempp dunecommon dunegrid)
target_link_libraries(traversal_benchmark bempp dunecommon dunegrid)
target_link_libraries(pure_traversal_benchmark bempp dunecommon dunegrid)

