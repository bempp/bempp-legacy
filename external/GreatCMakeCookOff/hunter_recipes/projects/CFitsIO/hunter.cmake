# Copyright (c) 2015, Damien Buhl 
# All rights reserved.

if(DEFINED HUNTER_CMAKE_PROJECTS_CFITSIO_HUNTER_CMAKE)
  return()
else()
  set(HUNTER_CMAKE_PROJECTS_CFITSIO_HUNTER_CMAKE 1)
endif()

include(hunter_add_version)
include(hunter_download)
include(hunter_pick_scheme)

hunter_add_version(
    PACKAGE_NAME
    CFitsIO
    VERSION
    "3350"
    URL
    "ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3350.tar.gz"
    SHA1
    e928832708d6a5df21a1e17ae4a63036cab7c1b9
)

hunter_add_version(
    PACKAGE_NAME
    CFitsIO
    VERSION
    "3360"
    URL
    "ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3360.tar.gz"
    SHA1
    0d3099721a6bf1e763f158c447a74c5e8192412d
)

hunter_add_version(
    PACKAGE_NAME
    CFitsIO
    VERSION
    "3370"
    URL
    "ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3370.tar.gz"
    SHA1
    48bd6389dcff3228508eec70384f2cae3a88ff32
)

hunter_pick_scheme(DEFAULT url_sha1_autotools)
hunter_download(PACKAGE_NAME CFitsIO)
