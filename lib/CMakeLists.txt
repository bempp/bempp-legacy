#Disable Numpy deprecated warning
add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
# Make CGAL Thread safe
add_definitions(-DCGAL_HAS_THREADS)

#Configure All Option files
foreach(config_file trilinos ahmed opencl data_types blas_and_lapack python)
    set(filename common/config_${config_file}.hpp)
    configure_file(${filename}.in
        ${PROJECT_BINARY_DIR}/include/bempp/${filename}
    )
endforeach()
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
    cuda_add_library(libbempp SHARED
        ${LIB_SOURCES} ${CUDA_SOURCES} ${LIB_HEADERS})
else ()
    add_library(libbempp SHARED ${LIB_SOURCES} ${LIB_HEADERS})
endif ()

# Target name already includes lib prefix
# However, it makes it easy to differentiate from python targets
set_target_properties(libbempp PROPERTIES PREFIX "")

# Add targets that may not yet exist, if they are being built as lookups
set(imported_targets
    # Dune targets
    dunecommon dunegeometry dunegrid 
)
unset(library_files)
foreach(library ${imported_targets})
    if(TARGET ${library})
        target_link_libraries(libbempp ${library})
    endif()
endforeach()

# Link with CGAL
target_link_libraries(libbempp ${CGAL_LIBRARY} ${CGAL_Core_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
    ${GMP_LIBRARIES} ${MPFR_LIBRARIES})



if (TBB_LIBRARY)
	target_link_libraries(libbempp optimized ${TBB_LIBRARY} ${TBB_MALLOC_LIBRARY})
endif()

# Link Cairo
# target_link_libraries(libbempp ${CAIRO_LIBRARIES})

if(TBB_LIBRARY_DEBUG)
    target_link_libraries(libbempp
        debug ${TBB_LIBRARY_DEBUG} ${TBB_MALLOC_LIBRARY_DEBUG})
endif()

if("${CMAKE_VERSION}" VERSION_LESS "2.8.11")
    set_target_properties(libbempp
        PROPERTIES COMPILE_DEFINITIONS_DEBUG TBB_USE_DEBUG=1)
else()
    target_compile_definitions(libbempp
        PRIVATE $<$<CONFIG:Debug>:TBB_USE_DEBUG=1>)
endif()
change_tbb_install_name(libbempp)

target_link_libraries(libbempp ${PYTHON_LIBRARY})
target_link_libraries(libbempp ${ZLIB_LIBRARIES})

if (WITH_CUDA)
   include(FindCUDA)
endif ()

# Install library
install(TARGETS libbempp
  EXPORT BemppTargets
  RUNTIME DESTINATION "${RUNTIME_INSTALL_PATH}"
  LIBRARY DESTINATION "${LIBRARY_INSTALL_PATH}"
  ARCHIVE DESTINATION "${LIBRARY_INSTALL_PATH}"
)

# Make sure its built after the external projects
depends_on_lookups(libbempp)

# Adds a target to copy headers
# Makes it possible to export the build directory,
# and have other projects build directly from the build tree.
include(TargetCopyFiles)
add_copy_directory(copy_headers
    "${CMAKE_CURRENT_SOURCE_DIR}"
    DESTINATION "${PROJECT_BINARY_DIR}/include/bempp"
    GLOB *.h *.hpp *.cl.str
    EXCLUDE *.in.h *.in.hpp
)
add_dependencies(libbempp copy_headers)

# Install headers
set(incpath "${INCLUDE_INSTALL_PATH}/bempp")
install(DIRECTORY "${PROJECT_SOURCE_DIR}/lib/"
    DESTINATION "${incpath}"
    FILES_MATCHING PATTERN "*.hpp"
)
install(FILES ${COMMON_HEADERS} DESTINATION "${incpath}/common")
install(FILES ${CL_FILES} DESTINATION "${incpath}/fiber/CL")
install(FILES ${CL_STRING_FILES} DESTINATION "${incpath}/fiber/CL")

# Add a target to generate API documentation with Doxygen
if (DOXYGEN_FOUND)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in"
        "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" @ONLY)
    file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/doc/cpp")
    add_custom_target(cpp_documentation
        ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/doc/cpp"
        COMMENT "Generate API documentation with Doxygen" VERBATIM
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
    )
    install(DIRECTORY "${PROJECT_BINARY_DIR}/doc/cpp"
        DESTINATION ${DOC_INSTALL_PATH}/bempp/cpp)
    add_dependencies(documentation cpp_documentation)
endif(DOXYGEN_FOUND)
