#Configure All Option files
configure_file(
	${CMAKE_SOURCE_DIR}/lib/common/config_trilinos.hpp.in
        ${CMAKE_BINARY_DIR}/include/bempp/common/config_trilinos.hpp)
configure_file(
	${CMAKE_SOURCE_DIR}/lib/common/config_ahmed.hpp.in
        ${CMAKE_BINARY_DIR}/include/bempp/common/config_ahmed.hpp)
configure_file(
	${CMAKE_SOURCE_DIR}/lib/common/config_opencl.hpp.in
        ${CMAKE_BINARY_DIR}/include/bempp/common/config_opencl.hpp)
configure_file(
        ${CMAKE_SOURCE_DIR}/lib/common/config_alugrid.hpp.in
        ${CMAKE_BINARY_DIR}/include/bempp/common/config_alugrid.hpp)
configure_file(
        ${CMAKE_SOURCE_DIR}/lib/common/config_data_types.hpp.in
        ${CMAKE_BINARY_DIR}/include/bempp/common/config_data_types.hpp)
# TODO: check that at most one of WITH_MKL, WITH_GOTOBLAS and WITH_OPENBLAS is set
configure_file(
        ${CMAKE_SOURCE_DIR}/lib/common/config_blas_and_lapack.hpp.in
        ${CMAKE_BINARY_DIR}/include/bempp/common/config_blas_and_lapack.hpp)
include_directories(${CMAKE_BINARY_DIR}/include)

# Find all source and header files
file(GLOB_RECURSE LIB_SOURCES *.cpp)
file(GLOB_RECURSE LIB_HEADERS *.hpp)
file(GLOB COMMON_HEADERS ${CMAKE_BINARY_DIR}/include/bempp/common/*.hpp)
file(GLOB CL_FILES ${CMAKE_SOURCE_DIR}/lib/fiber/CL/*.cl)
file(GLOB CL_STRING_FILES ${CMAKE_SOURCE_DIR}/lib/fiber/CL/*.str)

# Headers are ignored by make, but their list is used by IDEs
# like Qt Creator and Visual Studio to generate project file list.
if (WITH_CUDA)
   file(GLOB_RECURSE CUDA_SOURCES *.cu)
   cuda_add_library(bempp SHARED 
	${LIB_SOURCES}
	${CUDA_SOURCES} 
	${LIB_HEADERS}
	${DUNE_HEADERS}
	)
else ()
     add_library(bempp SHARED 
	${LIB_SOURCES} 
	${LIB_HEADERS}
	${DUNE_HEADERS}
	)
endif ()

# Add includes and target links

# BLAS and LAPACK
target_link_libraries(bempp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
include_directories(${BLAS_INCLUDE_DIR} ${LAPACK_INCLUDE_DIR})

# Boost
include_directories(${BOOST_INCLUDE_DIR})

# TBB
include_directories(${TBB_INCLUDE_DIR})
target_link_libraries (bempp debug ${TBB_LIBRARY_DEBUG} optimized ${TBB_LIBRARY})
set_property(
   TARGET bempp
   PROPERTY COMPILE_DEFINITIONS_DEBUG TBB_USE_DEBUG=1
   )

# CUDA
if (WITH_CUDA)
   include(FindCUDA)
endif ()

# AHMED
if (WITH_AHMED)
    target_link_libraries (bempp ${AHMED_LIB})
    include_directories(${AHMED_INCLUDE_DIR})
endif ()

# Dune
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/include)
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
  EXPORT BemppLibraryDepends
  RUNTIME DESTINATION bempp/bin
  LIBRARY DESTINATION bempp/lib
  ARCHIVE DESTINATION bempp/lib)

# Install CMake files for the library
install(EXPORT BemppLibraryDepends DESTINATION bempp/lib/cmake/Bempp)
configure_file(
	${CMAKE_SOURCE_DIR}/cmake/BemppConfig.cmake.in
        ${CMAKE_BINARY_DIR}/lib/cmake/Bempp/BemppConfig.cmake)


# Install headers
install(DIRECTORY ${CMAKE_SOURCE_DIR}/lib/ DESTINATION bempp/include/bempp FILES_MATCHING PATTERN "*.hpp")
install(FILES ${COMMON_HEADERS} DESTINATION bempp/include/bempp/common)
install(FILES ${CL_FILES} DESTINATION bempp/include/bempp/fiber/CL)
install(FILES ${CL_STRING_FILES} DESTINATION bempp/include/bempp/fiber/CL)
install(FILES ${CMAKE_BINARY_DIR}/lib/cmake/Bempp/BemppConfig.cmake DESTINATION bempp/lib/cmake/Bempp)


# Create correct rpath for installation
SET(BEMPP_RPATH "\$ORIGIN")

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

SET_TARGET_PROPERTIES(bempp
	PROPERTIES INSTALL_RPATH "${BEMPP_RPATH}")

# Do not explicitly link executables using libbempp with internal
# dependencies of BEM++ (like the Trilinos libraries, AHMED, BLAS...)
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
SET_TARGET_PROPERTIES(bempp
	PROPERTIES LINK_INTERFACE_LIBRARIES "")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")	
