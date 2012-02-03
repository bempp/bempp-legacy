# Find all source and header files
file(GLOB_RECURSE LIB_SOURCES *.cpp)
file(GLOB_RECURSE LIB_HEADERS *.hpp)
file(GLOB_RECURSE DUNE_FOAMGRID_HEADERS ${CMAKE_SOURCE_DIR}/contrib/dune/dune-foamgrid/dune/foamgrid/*.hh)
file(GLOB_RECURSE DUNE_GRID_HEADERS ${CMAKE_SOURCE_DIR}/contrib/dune/dune-grid/dune/grid/*.hh)
file(GLOB_RECURSE DUNE_COMMON_HEADERS ${CMAKE_SOURCE_DIR}/contrib/dune/dune-common/dune/common/*.hh)

add_definitions(-DCOMPILE_FOR_DOUBLE)
add_definitions(-DCHECK_ARRAY_BOUNDS)

# Headers are ignored by make, but their list is used by IDEs
# like Qt Creator and Visual Studio to generate project file list.
add_library(bempp SHARED 
	${LIB_SOURCES} 
	${LIB_HEADERS}
	${DUNE_FOAMGRID_HEADERS}
	${DUNE_GRID_HEADERS}
	${DUNE_COMMON_HEADERS}
	)

target_link_libraries (bempp
	${BLAS_LIB}
	${LIB_DUNE_COMMON}
	${LIB_DUNE_GRID}
	)


# Add a target to generate API documentation with Doxygen
find_package(Doxygen)
if (DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in 
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../doc)
  add_custom_target(doc 
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../doc
    COMMENT "Generate API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

# Install library
install(TARGETS bempp
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

# Install header files
install(FILES ${LIB_HEADERS}
  DESTINATION include/bempp)
