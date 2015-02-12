# - Try to find Umfpack
# Once done this will define
#
#  Umfpack_FOUND        - system has UMFPACK
#  Umfpack_INCLUDE_DIR  - include directories for UMFPACK
#  Umfpack_LIBRARIES    - libraries for UMFPACK

if (NOT Umfpack_LIBRARY_DIR)
    set(Umfpack_LIBRARY_DIR "${PROJECT_BINARY_DIR}/external/lib")
endif()

if (NOT Umfpack_INCLUDE_DIR)
    set(Umfpack_INCLUDE_DIR "${PROJECT_BINARY_DIR}/external/include")
endif()

# Check for header file
find_path(Umfpack_INCLUDE_DIR umfpack.h
 HINTS ${Umfpack_INCLUDE_DIR}
 DOC "Directory where the UMFPACK header is located"
 )

# Check for UMFPACK library
find_library(UMFPACK_LIBRARY umfpack
  HINTS ${Umfpack_LIBRARY_DIR}
  DOC "The UMFPACK library"
  )


# Check for AMD library
find_library(AMD_LIBRARY amd
  HINTS ${Umfpack_LIBRARY_DIR}
  DOC "The AMD library"
  )

# Check for CHOLMOD
find_library(CHOLMOD_LIBRARY cholmod
  HINTS ${Umfpack_LIBRARY_DIR}
  DOC "The CHOLMOD library"
)

# Check for COLAMD
find_library(COLAMD_LIBRARY colamd
  HINTS ${Umfpack_LIBRARY_DIR}
  DOC " The COLAMD Library"
)

# Check for CCOLAMD
find_library(CCOLAMD_LIBRARY ccolamd
  HINTS ${Umfpack_LIBRARY_DIR}
  DOC " The CCOLAMD Library"
)

# Check for CAMD
find_library(CAMD_LIBRARY camd
  PATHS ${Umfpack_LIBRARY_DIR}
  DOC " The CAMD Library"
)

# Collect libraries

if (UMFPACK_LIBRARY)
    set(Umfpack_LIBRARIES ${UMFPACK_LIBRARY})
endif()

if (AMD_LIBRARY)
  set(Umfpack_LIBRARIES ${Umfpack_LIBRARIES} ${AMD_LIBRARY})
endif()

if (CHOLMOD_LIBRARY)
  set(Umfpack_LIBRARIES ${Umfpack_LIBRARIES} ${CHOLMOD_LIBRARY})
endif()

if (COLAMD_LIBRARY)
  set(Umfpack_LIBRARIES ${Umfpack_LIBRARIES} ${COLAMD_LIBRARY})
endif()

if (CCOLAMD_LIBRARY)
  set(Umfpack_LIBRARIES ${Umfpack_LIBRARIES} ${CCOLAMD_LIBRARY})
endif()

if (CAMD_LIBRARY)
  set(Umfpack_LIBRARIES ${Umfpack_LIBRARIES} ${CAMD_LIBRARY})
endif()


# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Umfpack
  "UMFPACK could not be found. Be sure to set UMFPACK_DIR."
  Umfpack_LIBRARIES Umfpack_INCLUDE_DIR AMD_LIBRARY CHOLMOD_LIBRARY COLAMD_LIBRARY CCOLAMD_LIBRARY CAMD_LIBRARY)

if(NOT Umfpack_FOUND)
    unset(Umfpack_LIBRARIES CACHE)
    unset(UMFPACK_LIBRARY CACHE)
    unset(Umfpack_INCLUDE_DIR CACHE)
    unset(AMD_LIBRARY CACHE)
    unset(CHOLMOD_LIBRARY CACHE)
    unset(COLAMD_LIBRARY CACHE)
    unset(CCOLAMD_LIBRARY CACHE)
    unset(CAMD_LIBRARY CACHE)
endif()

