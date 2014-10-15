# Pulls together different dune components, some which are exported, some which are not.

# Look for common library first
unset(arguments)
if(Dune_FIND_REQUIRED)
    list(APPEND arguments REQUIRED)
endif()
unset(quietly)
if(Dune_FIND_QUIETLY)
    list(APPEND arguments QUIET)
    set(quietly QUIET)
endif()
find_package(dune-common ${arguments} ${quietly}
    HINTS ${dune_PREFIX} $ENV{dune_PREFIX})
if(NOT dune-common_FOUND)
    return()
endif()


# Remove foamgrid from components since provided within BEM
set(components ${Dune_FIND_COMPONENTS})
list(REMOVE_ITEM components foamgrid)

unset(required_not_found)
foreach(component ${components})
    find_package(dune-${component} ${quietly}
        HINTS ${dune_PREFIX} $ENV{dune_PREFIX})
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
