# Find all source and header files
file(GLOB_RECURSE LIB_SOURCES *.cpp)
file(GLOB_RECURSE LIB_HEADERS *.hpp)
file(GLOB_RECURSE DUNE_FOAMGRID_HEADERS ${CMAKE_SOURCE_DIR}/contrib/dune/dune-foamgrid/dune/foamgrid/*.hh)
file(GLOB_RECURSE DUNE_GRID_HEADERS ${CMAKE_SOURCE_DIR}/contrib/dune/dune-grid/dune/grid/*.hh)
file(GLOB_RECURSE DUNE_COMMON_HEADERS ${CMAKE_SOURCE_DIR}/contrib/dune/dune-common/dune/common/*.hh)
find_file(ARMA_HEADER armadillo ${ARMADILLO_INCLUDE_DIRS})

add_definitions(-DCOMPILE_FOR_DOUBLE)

# Headers are ignored by make, but their list is used by IDEs
# like Qt Creator and Visual Studio to generate project file list.
add_library(bempp SHARED 
	${LIB_SOURCES} 
	${LIB_HEADERS}
	${DUNE_FOAMGRID_HEADERS}
	${DUNE_GRID_HEADERS}
	${DUNE_COMMON_HEADERS}
	${ARMA_HEADER}
	)

target_link_libraries (bempp
	${BLAS_LIBRARIES}
	${LIB_DUNE_COMMON}
	${LIB_DUNE_GRID}
	)

# TBB

if (TBB_LIBRARY_DEBUG_FOUND)
	target_link_libraries (bempp debug ${TBB_LIBRARY_DEBUG} optimized ${TBB_LIBRARY_RELEASE})
else ()
	target_link_libraries (bempp ${TBB_LIBRARY_RELEASE})
endif ()

set_property(
   TARGET bempp
   PROPERTY COMPILE_DEFINITIONS_DEBUG TBB_USE_DEBUG=1
   )

include_directories(${TBB_INCLUDE_DIR})

if (WITH_AHMED)
    target_link_libraries (bempp ${AHMED_LIBRARY} ${METIS_LIBRARY})
    add_definitions(-DWITH_AHMED)
    include_directories(${AHMED_INCLUDE_DIR})
endif ()

if (WITH_TRILINOS)
    target_link_libraries (bempp ${EPETRA_LIBRARY})
    add_definitions(-DWITH_TRILINOS)
    include_directories(${TRILINOS_INCLUDE_DIR})
endif ()

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
