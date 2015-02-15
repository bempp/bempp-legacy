# Looks for VTK. If not found, download and install it.
if(VTK_ARGUMENTS)
    cmake_parse_arguments(VTK
        "URL;MD5;BUILD_TYPE;INSTALL_PREFIX"
        ${VTK_ARGUMENTS}
    )
endif()
if(NOT VTK_URL)
    set(arguments
        URL;
        http://www.vtk.org/files/release/6.1/VTK-6.1.0.tar.gz
        URL_MD5; 25e4dfb3bad778722dcaec80cd5dab7d 
    )
elseif(VTK_MD5)
    set(arguments URL;${VTK_URL};URL_MD5;${VTK_MD5})
else()
    message(FATAL_ERROR "URL specified, but no MD5. Aborting")
endif()
if(NOT VTK_BUILD_TYPE)
    set(VTK_BUILD_TYPE Release)
endif()

# Create a file of variables which is parsed by cmake when configuring VTK.
include(PassonVariables)
passon_variables(Trilinos
    FILENAME "${EXTERNAL_ROOT}/src/VtkVariables.cmake"
    PATTERNS
        "CMAKE_[^_]*_R?PATH" 
)
file(APPEND "${EXTERNAL_ROOT}/src/VtkVariables.cmake"
    "\nset(CMAKE_INSTALL_PREFIX \"${EXTERNAL_ROOT}\" CACHE PATH \"\")\n"
    "\nset(VTK_REQUIRED_OBJCXX_FLAGS \"\" CACHE STRING \"Extra flags for Objective C compilation\")\n"
)

set(VTK_CXX_COMPILER)
set(VTK_C_COMPILER)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(VTK_CXX_COMPILER /usr/bin/clang++)
    set(VTK_C_COMPILER /usr/bin/clang)
else()
    set(VTK_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
    set(VTK_C_COMPILER "${CMAKE_C_COMPILER}")
endif()


ExternalProject_Add(
    VTK 
    PREFIX ${EXTERNAL_ROOT}
    ${arguments}
    CMAKE_ARGS -DCMAKE_C_COMPILER=${VTK_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${VTK_CXX_COMPILER}
               -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
               -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
               -DBUILD_SHARED_LIBS:BOOL=ON
               -DCMAKE_BUILD_TYPE=${VTK_BUILD_TYPE}
               -DCMAKE_C_COMPILER=/usr/bin/clang
               -DCMAKE_CXX_COMPILER=/usr/bin/clang++
               -C ${EXTERNAL_ROOT}/src/VtkVariables.cmake
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

add_recursive_cmake_step(VTK DEPENDEES install)
