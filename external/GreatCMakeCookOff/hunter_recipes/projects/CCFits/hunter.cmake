# Copyright (c) 2015, Damien Buhl 
# All rights reserved.

if(DEFINED HUNTER_CMAKE_PROJECTS_CCFITS_HUNTER_CMAKE)
  return()
else()
  set(HUNTER_CMAKE_PROJECTS_CCFITS_HUNTER_CMAKE 1)
endif()

include(hunter_add_version)
include(hunter_download)
include(hunter_pick_scheme)
include(hunter_add_package)

# Makes it possible to use syste cfitsio
if(NOT CFitsIO_FOUND)
  hunter_add_package(CFitsIO)
  find_package(CFitsIO REQUIRED)
  if(NOT CFitsIO_FOUND)
    message(FATAL_ERROR "need cfitsio")
  endif()
endif()
get_filename_component(libcfitsio_dir "${CFitsIO_LIBRARY}" DIRECTORY)

hunter_add_version(
    PACKAGE_NAME
    CCFits
    VERSION
    "2.4"
    URL
    "http://heasarc.gsfc.nasa.gov/fitsio/ccfits/CCfits-2.4.tar.gz"
    SHA1
    3de2a6379bc1024300befae95cfdf33645a7b64a
)

hunter_pick_scheme(DEFAULT ccfits)
hunter_download(PACKAGE_NAME CCFits PACKAGE_DEPENDS_ON CFitsIO)
