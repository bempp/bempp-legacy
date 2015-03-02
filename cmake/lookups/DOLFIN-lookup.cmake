# Install Dolfin

if(NOT TARGET DOLFIN)
    include(PassonVariables)
    passon_variables(DOLFIN
        FILENAME "${EXTERNAL_ROOT}/src/DolfinVariables.cmake"
        PATTERNS
        "CMAKE_[^_]*_R?PATH" "CMAKE_C_.*" "CMAKE_CXX_.*"
        "BLAS_.*" "LAPACK_.*" "SWIG_*"
        )
endif()

if(NOT DOLFIN_BUILD_TYPE)
    set(DOLFIN_BUILD_TYPE Release)
endif()

unset(depends)
set(depends)
foreach(component Boost SWIG FFC Umfpack VTK Eigen3)
        if(TARGET ${component})
            list(APPEND depends ${component})
        endif()
endforeach()

set(DOLFIN_OMP)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(DOLFIN_OMP OFF)
else()
    set(DOLFIN_OMP ON)
endif()

include(PatchScript)
set(patchdir "${PROJECT_SOURCE_DIR}/cmake/patches/DOLFIN")
create_patch_script(DOLFIN patch_script
    CMDLINE "-p0"
    WORKING_DIRECTORY "${EXTERNAL_ROOT}/src/DOLFIN"
    "${patchdir}/dolfin_cmake.patch"
)

ExternalProject_Add(
    DOLFIN
    DEPENDS ${depends}
    URL http://bitbucket.org/fenics-project/dolfin/downloads/dolfin-1.5.0.tar.gz
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS ${depends}
    CONFIGURE_COMMAND /bin/bash -c "PYTHONPATH=${EXTERNAL_ROOT}/python cmake -D CMAKE_BUILD_TYPE=${DOLFIN_BUILD_TYPE} -D UFC_DIR=${EXTERNAL_ROOT}/share/ufc -D EIGEN3_INCLUDE_DIR=${Eigen3_INCLUDE_DIR} -D CMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT} -D DOLFIN_INSTALL_PYTHON_MODULE_DIR=${EXTERNAL_ROOT}/python -D DOLFIN_INSTALL_PYTHON_PURE_MODULE_DIR=${EXTERNAL_ROOT}/python -D DOLFIN_ENABLE_PETSC4PY:BOOL=OFF -D DOLFIN_SKIP_BUILD_TESTS:BOOL=ON -D DOLFIN_ENABLE_PETSC:BOOL=OFF -D DOLFIN_ENABLE_SLEPC4PY:BOOL=OFF -D DOLFIN_ENABLE_SLEPC:BOOL=OFF -D DOLFIN_ENABLE_OPENMP:BOOL=${DOLFIN_OMP} -D DOLFIN_ENABLE_MPI:BOOL=OFF -D PY_PLY_FOUND:BOOL=ON -C ${EXTERNAL_ROOT}/src/DolfinVariables.cmake ${EXTERNAL_ROOT}/src/DOLFIN/"
    PATCH_COMMAND ${patch_script}
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    )

add_recursive_cmake_step(DOLFIN DEPENDEES install)

