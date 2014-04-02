# Trie and finds dune
# Look for common library first
find_library(Dune_common_LIBRARY NAMES dunecommon DOC "Common dune functionalities")
if(NOT Dune_common_LIBRARY)
  message(STATUS "[dune] common library not found")
  return()
endif()
set(Dune_LIBRARIES ${Dune_common_LIBRARY})

# Look for common include second
find_path(Dune_INCLUDE_DIR array.hh PATH_SUFFIXES dune/common)
if(NOT Dune_INCLUDE_DIR)
  message(STATUS "[dune] common include path not found")
  return()
endif()

# Then look for extra components
# Get directories where to look for stuff
get_filename_component(Dune_LIBRARY_DIR ${Dune_common_LIBRARY} PATH)
get_filename_component(Dune_INCLUDE_COMPONENT_DIR ${Dune_INCLUDE_DIR} PATH)
set(LIBRARY_COMPONENTS Dune_common_LIBRARY)
# Some components have libraries, others not
set(LIBRARIZED_COMPONENTS geometry grid)
foreach(COMPONENT ${Dune_FIND_COMPONENTS})
  # First look for library, if it is expected to exist
  list(FIND LIBRARIZED_COMPONENTS ${COMPONENT} DO_LIBRARY)
  if(DO_LIBRARY GREATER -1)
    find_library(
      Dune_${COMPONENT}_LIBRARY NAMES dune${COMPONENT}
      PATHS ${Dune_LIBRARY_DIR}
      DOC "dune component library ${COMPONENT}."
      NO_DEFAULT_PATH
    )
    set(LIBRARY_COMPONENTS ${LIBRARY_COMPONENTS} Dune_${COMPONENT}_LIBRARY)
    set(Dune_${COMPONENT}_FOUND TRUE)
  endif()
  # Then looks for include directory
  unset(Dune_${COMPONENT}_FOUND)
  if(IS_DIRECTORY ${Dune_INCLUDE_COMPONENT_DIR}/${COMPONENT})
      set(Dune_${COMPONENT}_FOUND TRUE)
  endif()
endforeach()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Dune
    REQUIRED_VARS Dune_LIBRARIES Dune_INCLUDE_DIR
    HANDLE_COMPONENTS
)

if(NOT Dune_FOUND)
  foreach(COMPONENT ${LIBRARY_COMPONENTS} Dune_INCLUDE_DIR Dune_LIBRARIES)
    unset(${COMPONENT} CACHE)
  endforeach()
endif()
