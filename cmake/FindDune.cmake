# Pulls together different Dune components, some which are exported, some which are not.

unset(quietly)
if(Dune_FIND_QUIETLY)
    set(quietly QUIET)
endif()

# List of additional places to look for
set(hints ${Dune_PREFIX} $ENV{Dune_PREFIX} ${dune_PREFIX} $ENV{dune_PREFIX})
# Remove devel from components since handled separately
set(components ${Dune_FIND_COMPONENTS})
list(REMOVE_ITEM components devel)

foreach(component common ${components})
    find_package(dune-${component} ${quietly} HINTS ${hints})
    set(Dune_${component}_FOUND ${dune-${component}_FOUND})
endforeach()

# Add devel component
list(FIND Dune_FIND_COMPONENTS "devel" has_devel_component)
if(has_devel_component GREATER -1 AND dune-common_FOUND)
    find_program(DuneProject_PROGRAM duneproject
        PATHS "${EXTERNAL_ROOT}/src/dune-common/bin"
        HINTS ${hints} ${Dune_common_PREFIX}
        PATH_SUFFIXES bin
    )
    if(DuneProject_PROGRAM)
        set(Dune_devel_FOUND TRUE)
    endif()
endif()

include(FindPackageHandleStandardArgs)
set(dummy TRUE)
find_package_handle_standard_args(Dune REQUIRED_VARS dune-common_FOUND HANDLE_COMPONENTS)
