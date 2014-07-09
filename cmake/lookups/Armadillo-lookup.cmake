if(Armadillo_ARGUMENTS)
    cmake_parse_arguments(Armadillo
        ""
        "URL;MD5;TIMEOUT"
        ""
        ${Armadillo_ARGUMENTS}
    )
endif()
if(NOT Armadillo_URL)
    set(arguments
        URL;
        http://sourceforge.net/projects/arma/files/armadillo-4.300.8.tar.gz
        URL_HASH;
        SHA256=1d1a77ad7a74e8b4a74d7b71c5e3fe488e26283907b316618de0a8558c60173a
    )
elseif(Armadillo_MD5)
    set(arguments URL;${Armadillo_URL};URL_HASH;MD5=${Armadillo_MD5})
else()
    message(FATAL_ERROR "URL specified, but no MD5. Aborting")
endif()
if(NOT Armadillo_BUILD_TYPE)
    set(Armadillo_BUILD_TYPE Release)
endif()
if(NOT Armadillo_TIMEOUT)
    set(Armadillo_TIMEOUT 15)
endif()

# Create list of library paths.
# They will be added to CMAKE_LIBRARY_PATH
set(library_dirs ${CMAKE_LIBRARY_PATH})
foreach(library ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
   get_filename_component(directory "${library}" PATH)
   list(FIND library_dirs "${directory}" index)
   if(NOT index EQUAL -1)
       list(APPEND library_dirs "${directory}")
   endif()
endforeach()
string(REGEX REPLACE ";" " " library_dirs "${library_dirs}")
string(REGEX REPLACE ";" " " rpath_dirs "${CMAKE_INSTALL_RPATH}")

# Create script to pass variables to Armadillo
include(PassonVariables)
passon_variables(Armadillo
    FILENAME "${EXTERNAL_ROOT}/src/ArmadilloVariables.cmake"
    PATTERNS
        "CMAKE_[^_]*_R?PATH" "CMAKE_C_.*" "CMAKE_CXX_.*"
        "BLAS_.*" "LAPACK_.*"
    ALSOADD
        "\nset(CMAKE_INSTALL_PREFIX \"${EXTERNAL_ROOT}\" CACHE STRING \"\")\n"
        "set(CMAKE_LIBRARY_PATH ${library_dirs} CACHE PATH \"\" FORCE)\n"
        "set(CMAKE_INSTALL_RPATH ${rpath_dirs} CACHE INTERNAL \"\")\n"
)

include(PatchScript)
set(patchdir "${PROJECT_SOURCE_DIR}/cmake/patches/armadillo")
create_patch_script(armadillo patch_script
	    CMDLINE "-p0"
	    WORKING_DIRECTORY "${EXTERNAL_ROOT}/src/Armadillo"
	    "${patchdir}/armadillo_disable_hdf5.patch"
)


# Finally, we can create external project
ExternalProject_Add(
    Armadillo
    PREFIX ${EXTERNAL_ROOT}
    ${arguments}
    TIMEOUT ${Armadillo_TIMEOUT}
    CMAKE_ARGS
        -C "${EXTERNAL_ROOT}/src/ArmadilloVariables.cmake"
        -DCMAKE_BUILD_TYPE=Release
    TIMEOUT 10
    # Wrap download, configure and build steps in a script to log output
    PATCH_COMMAND ${patch_script}
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)
# Rerun cmake to capture new armadillo install
add_recursive_cmake_step(Armadillo DEPENDEES install)
