# Trie and finds dune

# Look for common library first
find_library(dune_common_LIBRARY NAMES dunecommon DOC "Common dune functionalities")
if(NOT dune_common_LIBRARY)
  message(STATUS "[dune] common library not found")
  return()
endif()
set(dune_LIBRARIES ${Dune_COMMON_LIBRARY})

# Look for common include second
find_path(dune_INCLUDE_DIR array.hh PATH_SUFFIXES dune/common)
if(NOT dune_INCLUDE_DIR)
  message(STATUS "[dune] common include path not found")
  return()
endif()

# Then look for extra components
# Get directories where to look for stuff
get_filename_component(dune_LIBRARY_DIR ${Dune_common_LIBRARY} PATH)
get_filename_component(dune_INCLUDE_COMPONENT_DIR ${Dune_INCLUDE_DIR} PATH)
set(LIBRARY_COMPONENTS dune_common_LIBRARY)
# Some components have libraries, others not
set(LIBRARIZED_COMPONENTS geometry grid)
foreach(COMPONENT ${dune_FIND_COMPONENTS})
  # First look for library, if it is expected to exist
  list(FIND LIBRARIZED_COMPONENTS ${COMPONENT} DO_LIBRARY)
  if(DO_LIBRARY GREATER -1)
    find_library( 
      dune_${COMPONENT}_LIBRARY NAMES dune${COMPONENT}
      PATHS ${dune_LIBRARY_DIR}
      DOC "dune component library ${COMPONENT}." 
      NO_DEFAULT_PATH
    )
    set(LIBRARY_COMPONENTS ${LIBRARY_COMPONENTS} dune_${COMPONENT}_LIBRARY)
    if("${dune_${COMPONENT}_LIBRARY}" STREQUAL "Dune_${COMPONENT}_LIBRARY-NOTFOUND")
      message(STATUS "[dune-${COMPONENT}] library NOT found")
    else()
      set(dune_LIBRARIES ${Dune_LIBRARIES} Dune_${COMPONENT}_LIBRARY)
      message(STATUS "[dune-${COMPONENT}] library found")
    endif()
  endif()
  # Then looks for include directory
  if(NOT IS_DIRECTORY ${dune_INCLUDE_COMPONENT_DIR}/${COMPONENT})
    unset(dune_INCLUDE_DIR CACHE)
    message(STATUS "[dune-${COMPONENT}] Include path not found")
  endif()
endforeach()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
  dune 
  DEFAULT_MSG 
  dune_LIBRARIES 
  dune_INCLUDE_DIR
  ${LIBRARY_COMPONENTS}
)

if(NOT dune_FOUND)
  foreach(COMPONENT ${LIBRARY_COMPONENTS} dune_INCLUDE_DIR Dune_LIBRARIES)
    unset(${COMPONENT} CACHE)
  endforeach()
endif()
