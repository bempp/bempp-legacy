# Copyright (c) 2015, Damien Buhl 
# All rights reserved.
include(hunter_add_version)
include(hunter_download)
include(hunter_cmake_args)
include(hunter_pick_scheme)
include(hunter_add_package)
include(hunter_configuration_types)

# Makes it possible to use system gmp
if(NOT GMP_FOUND)
  hunter_add_package(GMP)
  find_package(GMP REQUIRED)
  if(NOT GMP)
    message(FATAL_ERROR "need GMP")
  endif()
endif()

hunter_add_version(
    PACKAGE_NAME
    MPFR
    VERSION
    "3.1.4"
    URL
    "http://www.mpfr.org/mpfr-current/mpfr-3.1.4.tar.gz"
    SHA1
    272212c889d0ad6775ab6f315d668f3d01d1eaa3
)

hunter_pick_scheme(DEFAULT MPFR)
get_filename_component(GMP_LIB_DIR ${GMP_LIBRARIES} DIRECTORY)
hunter_cmake_args(MPFR CMAKE_ARGS GMP_INCLUDES="${GMP_INCLUDES}" GMP_LIB_DIR="${GMP_LIB_DIR}")
hunter_configuration_types(MPFR CONFIGURATION_TYPES Release)
hunter_download(PACKAGE_NAME MPFR PACKAGE_DEPENDS_ON GMP)
