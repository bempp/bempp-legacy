find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()
include(CachedVariables)

# Set some variables to delete later
set(testvar_one True CACHE BOOL "NOTHING")
set(testvar_two True CACHE BOOL "NOTHING")
set(a_testvar_one True CACHE BOOL "NOTHING")

# Deletes some of these variables
delete_variables("^testvar_.*")
if(testvar_one OR testvar_two)
    message(FATAL_ERROR "Could not delete variables")
endif()
if(NOT a_testvar_one)
    message(FATAL_ERROR "Deleted wrong variable")
endif()
