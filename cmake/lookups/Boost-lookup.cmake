set(toolset "")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(toolset "toolset=intel-linux")
    file(WRITE "${EXTERNAL_ROOT}/src/user-config.jam"
      "using intel : : \"${CMAKE_CXX_COMPILER}\" ; \n"
    )
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(toolset "toolset=gcc")
    file(WRITE "${EXTERNAL_ROOT}/src/user-config.jam"
      "using gcc : : \"${CMAKE_CXX_COMPILER}\" ; \n"
    )
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(toolset "toolset=clang")
    file(WRITE "${EXTERNAL_ROOT}/src/user-config.jam"
      "using clang : : \"${CMAKE_CXX_COMPILER}\" ; \n"
    )
else()
  message(FATAL_ERROR
     "Unknown compiler ${CMAKE_CXX_COMPILER_ID}."
     "Please install boost manually."
  )
endif()


file(WRITE "${PROJECT_BINARY_DIR}/CMakeFiles/external/boost_configure.sh"
    "#!${bash_EXECUTABLE}\n"
    "userjam=\"${EXTERNAL_ROOT}/src/user-config.jam\"\n"
    "[ -e $userjam ] && cp $userjam tools/build/src\n"
    "[ -e $userjam ] && cp $userjam tools/build/v2\n"
    "\n"
    "./b2 ${toolset} variant=release link=shared --with-test --with-filesystem --with-program_options --with-system --with-thread --with-iostreams\\\n"
    "    cxxflags=\"${CMAKE_CXX_FLAGS}\" --debug-configuration --verbose\n"
)
set(configure_command "${EXTERNAL_ROOT}/src/boost_configure.sh")
file(COPY "${PROJECT_BINARY_DIR}/CMakeFiles/external/boost_configure.sh"
    DESTINATION "${EXTERNAL_ROOT}/src"
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
)

find_program(PATCH_EXECUTABLE patch REQUIRED)
ExternalProject_Add(
    Boost
    PREFIX ${EXTERNAL_ROOT}
    # Downloads boost from url -- much faster than svn
    URL http://downloads.sourceforge.net/project/boost/boost/1.57.0/boost_1_57_0.tar.bz2
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ./bootstrap.sh
    BUILD_COMMAND ${configure_command}
    INSTALL_COMMAND ./b2 ${toolset} link=shared variant=release --with-test --with-filesystem --with-program_options --with-system --with-thread --with-iostreams
        --prefix=${EXTERNAL_ROOT} install
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)
# Rerun cmake to capture new boost install
add_recursive_cmake_step(Boost DEPENDEES install)
set(BOOST_ROOT "${EXTERNAL_ROOT}" CACHE INTERNAL "Prefix for Boost install")
# Makes sure those are not in the CACHE, otherwise, new version will not be found
unset(Boost_INCLUDE_DIR CACHE)
unset(Boost_LIBRARY_DIR CACHE)
