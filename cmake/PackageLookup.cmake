# Adds the ability to "look-up" a package
#
# This objective is to define a way to either find a package with find_package and/or,
# depending on choices and circumstances, download and install that package.
#
# The user should only have to include this file and add to their cmake files:
#
#     lookup_package(<name>)
#
# The name should match that of an existing `find_package(<name>)` file.
# The lookup_package function depends on files in directories in the cmake prefixes paths that have
# the name of the package:
#
# - ${CMAKE_MODULE_PATH}/${package}/${package}-lookup.cmake
# - ${CMAKE_MODULE_PATH}/${package}/LookUp${package}.cmake
# - ${CMAKE_MODULE_PATH}/${package}/lookup.cmake
# - ${CMAKE_MODULE_PATH}/${package}-lookup.cmake
# - ${CMAKE_MODULE_PATH}/LookUp${package}.cmake
# - ${CMAKE_MODULE_PATH}/${package}-lookup.cmake
# - ${CMAKE_LOOKUP_PATH}/${package}-lookup.cmake
# - ${CMAKE_LOOKUP_PATH}/LookUp${package}.cmake
#
# These files are included when the function lookup_package is called.
# The files will generally have the structure:
#
# ~~~cmake
# # Parses arguments specific to the lookup recipe
# # Optional step. Below, only a URL single-valued argument is specified.
# if(package_ARGUMENTS)
#     cmake_parse_arguments(package "" "URL" "" ${package_ARGUMENTS})
# else()
#     set(package_URL https://gaggledoo.doogaggle.com)
# endif()
# # The external project name `<package>` must match the package name exactly
# ExternalProject_Add(package
#   URL ${URL_
# )
# # Reincludes cmake so newly installed external can be found via find_package.
# # Optional step.
# add_recursive_cmake_step(...)
# ~~~
#
# This pattern will first attempt to find the package on the system. If it is not found, an external
# project to create it is added, with an extra step to rerun cmake and find the newly installed
# package.
#
# This file adds three functions/macro:
# 
# ~~~cmake
# lookup_package(<name>    # Name for find_package and lookup recipe files
#    [DOWNLOAD_BY_DEFAULT] # Always look it up as external project first.
#                          # This ensures the external project is always compiled specifically for
#                          # this project.
#    [ARGUMENTS <list>]    # Arguments specific to the look up recipe.
#                          # They will be available inside the recipe under the variable
#                          # ${name}_ARGUMENTS. Lookup recipes also have access to EXTERNAL_ROOT,
#                          # a variable specifying a standard location for external projects in the
#                          # build tree
#    [...]                 # Arguments passed on to `find_package`.
# )
# ~~~~
#
# ~~~~cmake
# # Adds a custom step to the external project that calls cmake recusively
# # This makes it possible for the newly built package to be installed.
# add_recursive_cmake_step(<name> # Still the same package name
#    <${name}_FOUND> # Variable set to true if package is found
#    [...]           # Passed on to ExternalProject_Add_Step
#                    # in general, it will be `DEPENDEES install`,
#                    # making this step the last.
# )
# ~~~
#
# ~~~~cmake
# # Makes sure TARGET is built after the looked up projects.
# depends_on_lookups(TARGET)
# ~~~
#

# Sets location where external project are included
if(NOT EXTERNAL_ROOT)
  set(EXTERNAL_ROOT ${CMAKE_BINARY_DIR}/external)
endif(NOT EXTERNAL_ROOT)
# Makes sure external projects are found by cmake
list(APPEND CMAKE_PREFIX_PATH ${EXTERNAL_ROOT})

# Adds a target for all external projects, so they can be made prior to others.
if(NOT TARGET lookup_dependencies)
    add_custom_target(lookup_dependencies ALL)
endif()

include(FindPackageHandleStandardArgs)
include(ExternalProject)

function(_find_lookup_recipe package OUTVAR)
    foreach(path ${CMAKE_MODULE_PATH})
      list(APPEND cmake_paths ${path}/${package})
    endforeach()
    set(LOOKUP_RECIPE LOOKUP_RECIPE-NOTFOUND)
    foreach(filename ${package}-lookup.cmake LookUp${package}.cmake)
        find_path(LOOKUP_RECIPE ${filename}
            PATHS ${CMAKE_LOOKUP_PATH} ${CMAKE_MODULE_PATH} ${cmake_paths}
            NO_DEFAULT_PATH
        )
        if(LOOKUP_RECIPE)
            set(${OUTVAR}_DIR "${LOOKUP_RECIPE}" PARENT_SCOPE)
            set(${OUTVAR}_FILE "${LOOKUP_RECIPE}/${filename}" PARENT_SCOPE)
            return()
        endif()
    endforeach()

    if(NOT LOOKUP_RECIPE)
        find_path(LOOKUP_RECIPE lookup.cmake
            PATHS ${cmake_paths}
            NO_DEFAULT_PATH
        )
    endif()

    if(LOOKUP_RECIPE)
      set(${OUTVAR}_DIR "${LOOKUP_RECIPE}" PARENT_SCOPE)
      set(${OUTVAR}_FILE "${LOOKUP_RECIPE}/lookup.cmake" PARENT_SCOPE)
    endif()
endfunction()

# Looks for a lookup package file and includes it.
macro(lookup_package package)
    set(solitos "DOWNLOAD_BY_DEFAULT;REQUIRED;QUIET;")
    set(multiplos "ARGUMENTS;COMPONENTS")
    cmake_parse_arguments(${package} "${solitos}" "" "${multiplos}" ${ARGN})

    # Reappends components
    if(${package}_COMPONENTS)
        list(APPEND ${package}_UNPARSED_ARGUMENTS COMPONENTS)
        list(APPEND ${package}_UNPARSED_ARGUMENTS ${${package}_COMPONENTS})
    endif()
    # Check whether recursive
    unset(recursive)
    set(REGEX REPLACE "(-|+| )" "_" SANENAME "${name}")
    if(${SANENAME}_RECURSIVE)
        set(recursive TRUE)
    endif()
    # First try and find package (unless downloading by default)
    set(dolook TRUE)
    if(${package}_DOWNLOAD_BY_DEFAULT AND NOT recursive)
        set(dolook FALSE)
    endif()
    # Figure out whether to add REQUIRED and QUIET keywords
    set(required "")
    if(recursive AND ${package}_REQUIRED)
        set(required REQUIRED)
    endif()
    set(quiet "")
    if(${package}_QIET)
        set(quiet QUIET)
    elseif(NOT recursive)
        set(quiet QUIET)
    endif()
    if(dolook)
        find_package(${package} ${${package}_UNPARSED_ARGUMENTS}
            ${required} ${quiet}
        )
    endif()
    # Sets lower and upper case versions.
    # Otherwise some package will be registered as not found.
    # This is a problem with changing cmake practices.
    string(TOUPPER "${package}" PACKAGE)
    if(${PACKAGE}_FOUND AND NOT "${package}" STREQUAL "${PACKAGE}")
        set(${package}_FOUND ${${PACKAGE}_FOUND})
    endif()
    # If package is not found, then look for a recipe to download and build it
    if(NOT ${package}_FOUND)
        _find_lookup_recipe(${package} ${package}_LOOKUP_RECIPE)
        if(NOT ${package}_LOOKUP_RECIPE_FILE)
            # Checks if package is required
            set(msg "Could not find recipe to lookup "
                    "${package} -- ${${package}_RECIPE_DIR}")
            if(${package}_REQUIRED)
                message(FATAL_ERROR ${msg})
            elseif(NOT ${package}_QUIET)
                message(STATUS ${msg})
            endif()
        else()
            if(NOT ${package}_QUIET AND NOT ${package}_DOWNLOAD_BY_DEFAULT)
                message(STATUS "Will attempt to download and install ${package}")
            elseif(NOT ${package}_QUIET)
                message(STATUS "Will download, build,"
                   " and install a local version of ${package}")
            endif()
            set(CURRENT_LOOKUP_DIRECTORY "${${package}_LOOKUP_RECIPE_DIR}")
            include(${${package}_LOOKUP_RECIPE_FILE})
            unset(CURRENT_LOOKUP_DIRECTORY)
            add_dependencies(lookup_dependencies ${package})
        endif()
    endif()
endmacro()

# Makes target depend on external dependencies
macro(depends_on_lookups TARGET)
    add_dependencies(${TARGET} lookup_dependencies)
endmacro()

# Adds an external step to an external project to rerun cmake
macro(add_recursive_cmake_step name check_var)
  set(REGEX REPLACE "(-|+| )" "_" SANENAME "${name}")
  set(cmake_arguments -DCMAKE_PROGRAM_PATH:PATH=${EXTERNAL_ROOT}/bin
                      -DCMAKE_LIBRARY_PATH:PATH=${EXTERNAL_ROOT}/lib
                      -DCMAKE_INCLUDE_PATH:PATH=${EXTERNAL_ROOT}/include
                      -DCMAKE_PREFIX_PATH:PATH=${EXTERNAL_ROOT}
                      -D${SANENAME}_RECURSIVE:BOOL=TRUE
                      --no-warn-unused-cli)
  if(NOT "${check_var}" STREQUAL "NOCHECK")
      set(cmake_arguments ${cmake_arguments} -D${SANENAME}_REQUIREDONRECURSE:INTERNAL=TRUE)
  endif()
  ExternalProject_Add_Step(
    ${name} reCMake
    COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR} ${cmake_arguments}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    ${ARGN}
  )
  if(${${SANENAME}_REQUIREDONRECURSE})
    if(NOT ${${check_var}})
      message(FATAL_ERROR "[${name}] Could not be downloaded and installed")
    else()
        message(STATUS "[${name}] ${check_var} ${${check_var}} Found")
    endif()
  endif()
  # Sets a variable saying we are building this source externally
  set(${name}_BUILT_AS_EXTERNAL_PROJECT TRUE)
endmacro()

# Avoids anoying cmake warning, by actually using the variables.
# The will be if the appropriate find_* is used. But won't be otherwise.
if(CMAKE_PROGRAM_PATH)
endif()
if(CMAKE_LIBRARY_PATH)
endif()
if(CMAKE_INCLUDE_PATH)
endif()

