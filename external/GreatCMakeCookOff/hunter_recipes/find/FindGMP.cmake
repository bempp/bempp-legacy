# Try to find the GNU Multiple Precision Arithmetic Library (GMP)
# See http://gmplib.org/

if (GMP_INCLUDES AND GMP_LIBRARIES)
  set(GMP_FIND_QUIETLY TRUE)
endif (GMP_INCLUDES AND GMP_LIBRARIES)

find_path(GMP_INCLUDES
  NAMES
  gmp.h
  PATHS
  $ENV{GMPDIR}
  ${INCLUDE_INSTALL_DIR}
)

find_library(GMP_LIBRARIES gmp PATHS $ENV{GMPDIR} ${LIB_INSTALL_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG
                                  GMP_INCLUDES GMP_LIBRARIES)
mark_as_advanced(GMP_INCLUDES GMP_LIBRARIES)

if(GMP_FOUND)
  if(GMP_LIBRARIES MATCHES "\\.a$")
    add_library(gmp STATIC IMPORTED GLOBAL)
  else()
    add_library(gmp SHARED IMPORTED GLOBAL)
  endif()
  set_target_properties(gmp PROPERTIES
      IMPORTED_LOCATION "${GMP_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDES}")
endif()
