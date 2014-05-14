#Configure All Option files
configure_file(
	${PROJECT_SOURCE_DIR}/lib/common/config_trilinos.hpp.in
        ${PROJECT_BINARY_DIR}/include/bempp/common/config_trilinos.hpp)
configure_file(
	${PROJECT_SOURCE_DIR}/lib/common/config_ahmed.hpp.in
        ${PROJECT_BINARY_DIR}/include/bempp/common/config_ahmed.hpp)
configure_file(
	${PROJECT_SOURCE_DIR}/lib/common/config_opencl.hpp.in
        ${PROJECT_BINARY_DIR}/include/bempp/common/config_opencl.hpp)
configure_file(
        ${PROJECT_SOURCE_DIR}/lib/common/config_alugrid.hpp.in
        ${PROJECT_BINARY_DIR}/include/bempp/common/config_alugrid.hpp)
configure_file(
        ${PROJECT_SOURCE_DIR}/lib/common/config_data_types.hpp.in
        ${PROJECT_BINARY_DIR}/include/bempp/common/config_data_types.hpp)
# TODO: check that at most one of WITH_MKL, WITH_GOTOBLAS and WITH_OPENBLAS is set
configure_file(
        ${PROJECT_SOURCE_DIR}/lib/common/config_blas_and_lapack.hpp.in
        ${PROJECT_BINARY_DIR}/include/bempp/common/config_blas_and_lapack.hpp)

configure_file(${PROJECT_SOURCE_DIR}/lib/config_bempp.in.hpp
    ${PROJECT_BINARY_DIR}/include/bempp/config_bempp.hpp)

# Find all source and header files
file(GLOB_RECURSE LIB_SOURCES *.cpp)
file(GLOB_RECURSE LIB_HEADERS *.hpp)
file(GLOB COMMON_HEADERS ${PROJECT_BINARY_DIR}/include/bempp/common/*.hpp)
file(GLOB CL_FILES ${PROJECT_SOURCE_DIR}/lib/fiber/CL/*.cl)
file(GLOB CL_STRING_FILES ${PROJECT_SOURCE_DIR}/lib/fiber/CL/*.str)

# Headers are ignored by make, but their list is used by IDEs
# like Qt Creator and Visual Studio to generate project file list.
if (WITH_CUDA)
   file(GLOB_RECURSE CUDA_SOURCES *.cu)
   cuda_add_library(bempp SHARED
        ${LIB_SOURCES}
        ${CUDA_SOURCES}
        ${LIB_HEADERS}
	)
else ()
    add_library(bempp SHARED ${LIB_SOURCES} ${LIB_HEADERS})
endif ()

# Add includes and target links
target_link_libraries (bempp
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${dune-common_STATIC_LDFLAGS}
    ${dune-grid_STATIC_LDFLAGS}
    armadillo

    # Trilinos libraries, imported by cmake as targets
    amesos stratimikos stratimikosbelos
    belos epetra thyracore
    teuchoscomm teuchoscore
    teuchosnumerics teuchosparameterlist
    teuchosremainder rtop
)
if(TBB_LIBRARY_DEBUG)
    target_link_libraries (bempp
        debug ${TBB_LIBRARY_DEBUG} ${TBB_MALLOC_LIBRARY_DEBUG}
        optimized ${TBB_LIBRARY} ${TBB_MALLOC_LIBRARY}
    )
endif()

set_property(
   TARGET bempp
   PROPERTY COMPILE_DEFINITIONS_DEBUG TBB_USE_DEBUG=1
)
if(ALUGrid_GITTER_INCLUDE_DIR)
    get_source_file_property(file_flags
        grid/grid_factory.cpp COMPILE_FLAGS
    )
    set_source_files_properties(
        grid/grid_factory.cpp
        PROPERTIES
        COMPILE_FLAGS "${file_flagS} -I${ALUGrid_GITTER_INCLUDE_DIR}"
    )
endif()

if (WITH_CUDA)
   include(FindCUDA)
endif ()

if (WITH_AHMED)
    target_link_libraries (bempp ${AHMED_LIB})
    include_directories(${AHMED_INCLUDE_DIR})
endif ()


# Install library
install(TARGETS bempp
  EXPORT BemppLibraryDepends
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

# Make sure its built after the external projects
depends_on_lookups(bempp)

# Install headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/lib/ DESTINATION include/bempp FILES_MATCHING PATTERN "*.hpp")
install(FILES ${COMMON_HEADERS} DESTINATION include/bempp/common)
install(FILES ${CL_FILES} DESTINATION include/bempp/fiber/CL)
install(FILES ${CL_STRING_FILES} DESTINATION include/bempp/fiber/CL)

# Add a target to generate API documentation with Doxygen
if (DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../doc)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../doc
    COMMENT "Generate API documentation with Doxygen" VERBATIM)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../doc
          DESTINATION doc/bempp)
endif(DOXYGEN_FOUND)
