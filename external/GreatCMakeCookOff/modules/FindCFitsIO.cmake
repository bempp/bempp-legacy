# Defines the following variables
#
# - CFitsIO_FOUND if the library is found
# - CFitsIO_LIBRARY is the path to the library
# - CFitsIO_INCLUDE_DIR is the path to the include directory
# - CFitsIO_VERSION_STRING is the version of the library
if(CFitsIO_FOUND)
  return()
endif()

find_library(CFitsIO_LIBRARY cfitsio DOC "Path to the cfitsio library")
if(NOT "$ENV{CASAPATH}" STREQUAL "")
    string(FIND "$ENV{CASAPATH}" " " endpath)
    string(SUBSTRING "$ENV{CASAPATH}" 0 ${endpath} casapath)

    find_library(
      CFitsIO_LIBRARY cfitsio
      NAMES libcfitsio${CMAKE_SHARED_LIBRARY_SUFFIX}.1 libcfitsio${CMAKE_SHARED_LIBRARY_SUFFIX}.0
      PATHS "${casapath}" "${casapath}/Frameworks"
      DOC "Path to the cfitsio library"
    )
endif()
find_path(
  CFitsIO_INCLUDE_DIR fitsio.h
  PATH_SUFFIXES include include/cfitsio Frameworks/Headers/cfitsio
  DOC "Path to the cfitsio include directory"
  PATHS "${casapath}"
)

if(CFitsIO_INCLUDE_DIR)
  file(
    STRINGS ${CFitsIO_INCLUDE_DIR}/fitsio.h
    CFitsIO_VERSION_STRING
    REGEX "#define[ ]+CFITSIO_VERSION[ ]+([0-9]*\\.[0-9]*)"
  )
  string(
    REGEX REPLACE
    ".*#define[ ]+CFITSIO_VERSION[ ]+([0-9]*\\.[0-9]*).*"
    "\\1"
    CFitsIO_VERSION_STRING
    "${CFitsIO_VERSION_STRING}"
  )
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
  CFitsIO
  REQUIRED_VARS CFitsIO_LIBRARY CFitsIO_INCLUDE_DIR
  VERSION_VAR CFitsIO_VERSION_STRING
)
if(CFITSIO_FOUND AND NOT CFitsIO_FOUND)
    set(CFitsIO_FOUND ${CFITSIO_FOUND})
endif()
if(CFitsIO_FOUND)
  if(CFitsIO_LIBRARY MATCHES "\\.a$")
    add_library(cfitsio STATIC IMPORTED GLOBAL)
  else()
    add_library(cfitsio SHARED IMPORTED GLOBAL)
  endif()
  set_target_properties(cfitsio PROPERTIES
    IMPORTED_LOCATION "${CFitsIO_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CFitsIO_INCLUDE_DIR}")
endif()
