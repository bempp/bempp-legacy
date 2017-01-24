find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()
include(CachedVariables)

if(DEFINED testvar_one)
  unset(testvar_one CACHE)
  unset(testvar_one)
endif()
if(DEFINED testvar_two)
  unset(testvar_two CACHE)
  unset(testvar_two)
endif()

cached_variables(nothing "testvar_.*")
if(NOT "${nothing}" STREQUAL "")
  message(FATAL_ERROR "Found an unexpected variable")
endif()

set(testvar_one True CACHE BOOL "NOTHING")
set(testvar_two True CACHE BOOL "NOTHING")

cached_variables(one_only "testvar_o.*")
if(NOT "${one_only}" STREQUAL "-Dtestvar_one=\"True\"")
    message(FATAL_ERROR "Could not find variable")
endif()

cached_variables(two_vars "testvar_.*")
list(FIND two_vars "-Dtestvar_one=\"True\"" one_present)
list(FIND two_vars "-Dtestvar_two=\"True\"" two_present)
if(one_present EQUAL -1 OR two_present EQUAL -1)
    message(FATAL_ERROR "Could not find both variable")
endif()

