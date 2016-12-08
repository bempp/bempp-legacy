if(Eigen3_ARGUMENTS)
    cmake_parse_arguments(Eigen3
        ""
        "HG_REPOSITORY;HG_TAG;URL;MD5;TIMEOUT"
        ""
        ${Eigen3_ARGUMENTS}
    )
endif()

if(NOT Eigen3_HG_REPOSITORY AND NOT Eigen3_URL)
  # Can't trust cmake to download file. It will fail if not compiled
  # with SSL
  set(file_found FALSE)
  set(Eigen3_URL "${EXTERNAL_ROOT}/eigen.tgz")
  set(Eigen3_MD5 6a578dba42d1c578d531ab5b6fa3f741)
  if(EXISTS "${EXTERNAL_ROOT}/eigen.tgz")
      file(MD5 "${EXTERNAL_ROOT}/eigen.tgz" file_md5)
      if(file_md5 EQUAL "${Eigen3_MD5}")
          set(file_found TRUE)
      endif()
  endif()
  if(NOT file_found)
      file(MAKE_DIRECTORY "${EXTERNAL_ROOT}")
      find_package(Wget)
      if(WGET_FOUND)
          execute_process(COMMAND ${WGET_EXECUTABLE}
            http://bitbucket.org/eigen/eigen/get/3.2.9.tar.gz
            -O ${EXTERNAL_ROOT}/eigen.tgz
          )
      else()
          find_program(CURL_EXECUTABLE curl)
          execute_process(COMMAND ${CURL_EXECUTABLE}
            -L http://bitbucket.org/eigen/eigen/get/3.2.9.tar.gz
            -o ${EXTERNAL_ROOT}/eigen.tgz
          )
      endif()
  endif()
endif()
if(NOT Eigen3_BUILD_TYPE)
    set(Eigen3_BUILD_TYPE Release)
endif()
if(NOT Eigen3_TIMEOUT)
    set(Eigen3_TIMEOUT 10)
endif()

if(NOT "${Eigen3_URL}" STREQUAL "")
    if("${Eigen3_MD5}" STREQUAL "")
        message(FATAL_ERROR "Downloading from an URL requires an MD5 hash")
    endif()
    set(arguments URL "${Eigen3_URL}" URL_MD5 ${Eigen3_MD5})
elseif("${CMAKE_VERSION}" VERSION_GREATER 2.8.9 )
    set(arguments HG_REPOSITORY "${Eigen3_HG_REPOSITORY}")
    if(NOT Eigen3_HG_TAG)
        list(APPEND arguments HG_TAG ${Eigen3_HG_TAG})
    endif()
else()
    message(FATAL_ERROR "Requires cmake >= 2.8.10 to download Eigen from mercurial")
endif()
list(APPEND arguments TIMEOUT ${Eigen3_TIMEOUT})

# write subset of variables to cache for eigen3 to use
include(PassonVariables)
passon_variables(Eigen3
    FILENAME "${EXTERNAL_ROOT}/src/Eigen3Variables.cmake"
    PATTERNS
        "CMAKE_[^_]*_R?PATH"
        "CMAKE_Fortran_.*"
        "CMAKE_C_.*"
        "CMAKE_CXX_.*"
        "BLAS_.*" "FFTW3_.*"
        ".*_INCLUDE_DIRS?"
    ALSOADD
        "\nset(CMAKE_INSTALL_PREFIX \"${EXTERNAL_ROOT}\" CACHE STRING \"\")\n"
)

# Finally add project
ExternalProject_Add(
  Lookup-Eigen3
    PREFIX ${EXTERNAL_ROOT}
    ${arguments}
    CMAKE_ARGS
        -C "${EXTERNAL_ROOT}/src/Eigen3Variables.cmake"
        -DBUILD_SHARED_LIBS=OFF
        -DCMAKE_BUILD_TYPE=${Eigen3_BUILD_TYPE}
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

add_recursive_cmake_step(Lookup-Eigen3 DEPENDEES install)
