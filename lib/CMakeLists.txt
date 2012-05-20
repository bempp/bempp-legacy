#Configure All Option files
configure_file(
	${CMAKE_SOURCE_DIR}/lib/common/config_trilinos.hpp.in
        ${CMAKE_BINARY_DIR}/include/common/config_trilinos.hpp)
configure_file(
	${CMAKE_SOURCE_DIR}/lib/common/config_ahmed.hpp.in
        ${CMAKE_BINARY_DIR}/include/common/config_ahmed.hpp)
configure_file(
	${CMAKE_SOURCE_DIR}/lib/common/config_opencl.hpp.in
        ${CMAKE_BINARY_DIR}/include/common/config_opencl.hpp)
configure_file(
        ${CMAKE_SOURCE_DIR}/lib/common/config_alugrid.hpp.in
        ${CMAKE_BINARY_DIR}/include/common/config_alugrid.hpp)
configure_file(
        ${CMAKE_SOURCE_DIR}/lib/common/config_data_types.hpp.in
        ${CMAKE_BINARY_DIR}/include/common/config_data_types.hpp)
include_directories(${CMAKE_BINARY_DIR}/include/common)

# Find all source and header files
file(GLOB_RECURSE LIB_SOURCES *.cpp)
file(GLOB_RECURSE LIB_HEADERS *.hpp)


# Headers are ignored by make, but their list is used by IDEs
# like Qt Creator and Visual Studio to generate project file list.
add_library(bempp SHARED 
	${LIB_SOURCES} 
	${LIB_HEADERS}
	${DUNE_FOAMGRID_HEADERS}
	${DUNE_GRID_HEADERS}
	${DUNE_COMMON_HEADERS}
        ${DUNE_LOCALFUNCTIONS_HEADERS}
	)

# Add includes and target links

# BLAS and LAPACK
target_link_libraries(bempp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Boost
include_directories(${BOOST_INCLUDE_DIR})

# TBB
include_directories(${TBB_INCLUDE_DIR})
target_link_libraries (bempp debug ${TBB_LIBRARY_DEBUG} optimized ${TBB_LIBRARY})
set_property(
   TARGET bempp
   PROPERTY COMPILE_DEFINITIONS_DEBUG TBB_USE_DEBUG=1
   )

# AHMED
if (WITH_AHMED)
    target_link_libraries (bempp ${AHMED_LIB} ${METIS_LIB})
    include_directories(${AHMED_INCLUDE_DIR})
endif ()

# Dune
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/dune-common)
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/dune-grid)
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/dune-localfunctions)
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/dune-foamgrid)
target_link_libraries (bempp
	${LIB_DUNE_COMMON}
	${LIB_DUNE_GRID}
	)

#Trilinos
include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
target_link_libraries (bempp amesos
                             stratimikos
                             stratimikosbelos
                             belos
                             epetra
                             thyracore
			     teuchos
			     rtop
)

# ARMADILLO
include_directories(${ARMADILLO_INCLUDE_DIR})

# Install library
install(TARGETS bempp
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

# Install header files
install(FILES ${LIB_HEADERS}
  DESTINATION include)
