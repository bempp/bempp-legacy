include_directories(${CMAKE_BINARY_DIR}/include)


function(find_subdirectories location outvar)
  file(GLOB ALL_FILES RELATIVE ${location} ${location}/*)
  set(RESULT_VAR) 
  foreach(candidate ${ALL_FILES})
    if(IS_DIRECTORY ${location}/${candidate})
      set(RESULT_VAR ${RESULT_VAR} ${candidate})
    endif()
  endforeach()
  set(${outvar} ${RESULT_VAR} PARENT_SCOPE)
endfunction(find_subdirectories)

add_definitions(-DBOOST_TEST_DYN_LINK)
#find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/include)

# Figure out environment variable with wich to run test
get_filename_component(ENVVARS ${BOOST_UNIT_TEST_LIB} PATH)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(ENVVARS "DYLD_LIBRARY_PATH=${ENVVARS}:$ENV{DYLD_LIBRARY_PATH}")
else()
  set(ENVVARS "LD_LIBRARY_PATH=${ENVVARS}:$ENV{LD_LIBRARY_PATH}")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Files to be included in all tests.
file(GLOB TEST_MAIN_SOURCES *.cpp)
file(GLOB TEST_MAIN_HEADERS *.hpp)
find_subdirectories(${CMAKE_CURRENT_SOURCE_DIR} SUBDIRECTORIES)
# meshes has no test
list(REMOVE_ITEM SUBDIRECTORIES "meshes")

foreach(directory ${SUBDIRECTORIES}) 

  # Files for this particular test
  file(GLOB TEST_SOURCES ${directory}/*.cpp)
  file(GLOB TEST_HEADERS ${directory}/*.hpp)
  if(${directory} STREQUAL "space")
    set(TEST_SOURCES ${TEST_SOURCES} assembly/create_regular_grid.cpp)
  endif()
  add_executable(test_${directory} ${TEST_MAIN_SOURCES} ${TEST_SOURCES})
    
  target_link_libraries(test_${directory} bempp ${BOOST_UNIT_TEST_LIB})

  add_test(NAME test_${directory} 
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND test_${directory} )
  set_tests_properties(test_${directory} PROPERTIES ENVIRONMENT ${ENVVARS})
endforeach()

# Meshes
file(GLOB_RECURSE TEST_MESHES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    *.msh)
if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    foreach (f ${TEST_MESHES})
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${f}
                       ${CMAKE_CURRENT_BINARY_DIR}/${f}
                       COPYONLY)
    endforeach ()
endif ()
