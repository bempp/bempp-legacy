# Create a script that cmake can call
# This should remove some issues that arises when cmake tries to build
# a complicated command line
file(WRITE "${PROJECT_BINARY_DIR}/CMakeFiles/external/alugrid_configure.sh"
    "#!${bash_EXECUTABLE}\n"
    "# Calls configure script for ALUGrid packages\n"
    "export CC=${CMAKE_C_COMPILER}\n"
    "export CFLAGS=\"${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}\"\n"
    "export CXX=${CMAKE_CXX_COMPILER}\n"
    "export CXXFLAGS=\"${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}\"\n"
    "export FC=${CMAKE_Fortran_COMPILER}\n"
    "export F77=${CMAKE_Fortran_COMPILER}\n"
    "export F90=${CMAKE_Fortran_COMPILER}\n"
    "export FCFLAGS=\"${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}\"\n"
    "export PKG_CONFIG_PATH=\"$ENV{PKG_CONFIG_PATH}\"\n"
    "\n"
    "./configure"
       " --disable-shared"
       " --enable-static"
       " --with-pic"
       " --prefix=\"${EXTERNAL_ROOT}\""
       " --includedir=\"${EXTERNAL_ROOT}/include/ALUGrid\""
)
set(configure_command "${EXTERNAL_ROOT}/src/alugrid_configure.sh")
file(COPY "${PROJECT_BINARY_DIR}/CMakeFiles/external/alugrid_configure.sh"
    DESTINATION "${EXTERNAL_ROOT}/src/"
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
)


ExternalProject_Add(
    ALUGrid
    PREFIX ${EXTERNAL_ROOT}
    DEPENDS ${depends_on}
    URL http://aam.mathematik.uni-freiburg.de/IAM/Research/alugrid/ALUGrid-1.52.tar.gz
    URL_HASH MD5=393c3d8ac1e9def4280765b16b95f5f1
    CONFIGURE_COMMAND ${configure_command}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

