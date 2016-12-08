# Adds subdirectory to CMAKE_MODULE_PATH so the scripts of the GreatCMakeCookOff can be found.
set(GREAT_CMAKE_COOKOFF_MODULE_DIR ${CMAKE_CURRENT_LIST_DIR}/../modules CACHE DOC
    "Path to GreatCMakeCookOff module directory")
set(GREAT_CMAKE_COOKOFF_SCRIPT_DIR ${CMAKE_CURRENT_LIST_DIR}/../scripts CACHE DOC
    "Path to GreatCMakeCookOff script directory")
set(GREAT_CMAKE_COOKOFF_HUNTER_RECIPES_DIR ${CMAKE_CURRENT_LIST_DIR}/../hunter_recipes CACHE DOC
    "Path to GreatCMakeCookOff script directory")
macro(initialize_cookoff)
    if(CMAKE_VERSION VERSION_LESS 2.8)
        message(FATAL_ERROR "The GreatCMakeCookOff requires CMake > 2.8")
    endif()
    cmake_policy(GET CMP0012 policy_twelve)
    if(NOT ${policy_twelve} STREQUAL "NEW")
        message(FATAL_ERROR "The GreatCMakeCookOff requires policy CMP0012 NEW")
    endif()
    list(APPEND CMAKE_MODULE_PATH ${GREAT_CMAKE_COOKOFF_MODULE_DIR}
        ${GREAT_CMAKE_COOKOFF_SCRIPT_DIR})
endmacro()
