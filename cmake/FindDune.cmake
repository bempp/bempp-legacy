# Pulls together different dune components, some which are exported, some which are not.

# Look for common library first
if(Dune_FIND_REQUIRED)
    list(APPEND arguments REQUIRED)
endif()
unset(quietly)
if(Dune_FIND_QUIETLY)
    list(APPEND arguments QUIET)
    set(quietly QUIET)
endif()
pkg_check_modules(dune-common ${arguments} dune-common)
if(NOT dune-common_FOUND)
    return()
endif()


get_filename_component(Dune_INCLUDE_COMPONENT_DIR
    ${dune-common_INCLUDE_DIRS} PATH)
unset(required_not_found)
foreach(component ${Dune_FIND_COMPONENTS})
    pkg_check_modules(dune-${component} ${quietly} dune-${component})
    if(Dune_FIND_REQUIRED_${component} AND NOT dune-${component}_FOUND)
        list(APPEND required_not_found ${component})
    endif()
endforeach()
if("${required_not_found}" STREQUAL "")
    set(Dune_FOUND TRUE)
elseif(Dune_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find dune components ${required_not_found}")
else()
    message(STATUS "Could not find dune components ${required_not_found}")
endif()
