# simple CMake file to install tbb from its extracted download file to places we want.
cmake_minimum_required(VERSION 2.8)

# Find patch files and apply
file(GLOB PATCH_FILES ${ROOT}/install/tbb_*.patch)
find_program(PATCH_EXECUTABLE patch)
foreach(patch ${PATCH_FILES})
  message(STATUS "[PATCH] ${PATCH_EXECUTABLE} ${patch}")
  execute_process(
    COMMAND ${PATCH_EXECUTABLE} ${patch}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
endforeach()

# find libraries
file(GLOB LIBRARIES lib/lib*.*)
# find headers
file(GLOB_RECURSE INCLUDES include/*.h)

# install stuff
install(FILES ${LIBRARIES} DESTINATION lib)
install(FILES ${INCLUDES} DESTINATION include)
