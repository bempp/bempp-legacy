# Exports BEMPP so other packages can access it
export(TARGETS libbempp 
    FILE "${PROJECT_BINARY_DIR}/BemppTargets.cmake")

# Avoids creating an entry in the cmake registry.
if(NOT NOEXPORT)
    export(PACKAGE Bempp)
endif()

# First in binary dir
set(ALL_INCLUDE_DIRS ${BEMPP_INCLUDE_DIRS})
set(BEMPP_CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
get_filename_component(BEMPP_PREFIX_PATH "${EXTERNAL_ROOT}" ABSOLUTE)
configure_File(cmake/BemppConfig.in.cmake
    "${PROJECT_BINARY_DIR}/BemppConfig.cmake" @ONLY
)
configure_File(cmake/BemppConfigVersion.in.cmake
    "${PROJECT_BINARY_DIR}/BemppConfigVersion.cmake" @ONLY
)

# Then for installation tree
set(share_path "${SHARE_INSTALL_PATH}/cmake/Bempp")
if(NOT IS_ABSOLUTE "${share_path}")
    set(share_path "${CMAKE_INSTALL_PREFIX}/${share_path}")
endif()
set(BEMPP_CMAKE_MODULE_PATH "${share_path}/cmake")
set(include_path "${INCLUDE_INSTALL_PATH}/Bempp")
if(NOT IS_ABSOLUTE "${include_path}")
    set(include_path "${CMAKE_INSTALL_PREFIX}/${include_path}")
endif()
file(RELATIVE_PATH REL_INCLUDE_DIR "${share_path}" "${include_path}")
set(ALL_INCLUDE_DIRS "\${Bempp_CMAKE_DIR}/${REL_INCLUDE_DIR}")
unset(BEMPP_PREFIX_PATH)
configure_File(cmake/BemppConfig.in.cmake
    "${PROJECT_BINARY_DIR}/CMakeFiles/BemppConfig.cmake" @ONLY
)

# Finally install all files
install(FILES
    "${PROJECT_BINARY_DIR}/CMakeFiles/BemppConfig.cmake"
    "${PROJECT_BINARY_DIR}/BemppConfigVersion.cmake"
    DESTINATION ${SHARE_INSTALL_PATH}/cmake/Bempp
    COMPONENT dev
)

install(EXPORT BemppTargets
    DESTINATION ${SHARE_INSTALL_PATH}/cmake/Bempp
    COMPONENT dev
)

