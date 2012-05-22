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
file(GLOB COMMON_HEADERS ${CMAKE_BINARY_DIR}/include/common/*.hpp)

# Headers are ignored by make, but their list is used by IDEs
# like Qt Creator and Visual Studio to generate project file list.
add_library(bempp SHARED 
	${LIB_SOURCES} 
	${LIB_HEADERS}
	${DUNE_HEADERS}
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
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/include)
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
  RUNTIME DESTINATION bempp/bin
  LIBRARY DESTINATION bempp/lib
  ARCHIVE DESTINATION bempp/lib)

# Install headers
install(DIRECTORY ${CMAKE_SOURCE_DIR}/lib/ DESTINATION bempp/include FILES_MATCHING PATTERN "*.hpp")
install(FILES ${COMMON_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/bempp/include/common)


# Create correct rpath for installation
SET(BEMPP_RPATH "\$ORIGIN/../contrib/dune/lib;\$ORIGIN/../contrib/trilinos/lib;\$ORIGIN/../contrib/boost/lib/boost-1.49.0")
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   # Add Tbb for Mac
   SET(BEMPP_RPATH "${BEMPP_RPATH};\$ORIGIN/../contrib/tbb/lib")
endif()
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   # Add Tbb for Linux
   SET(BEMPP_RPATH "${BEMPP_RPATH};\$ORIGIN/../contrib/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21")
endif()
if (WITH_AHMED)
   get_filename_component(AHMED_PATH "${AHMED_LIB}" PATH)
   get_filename_component(METIS_PATH "${METIS_LIB}" PATH)
   SET(BEMPP_RPATH "${BEMPP_RPATH};${AHMED_PATH};${METIS_PATH}")
endif()

SET_TARGET_PROPERTIES(bempp
	PROPERTIES INSTALL_RPATH "${BEMPP_RPATH}")
	
