# Adds the ability to "look-up" a package, e.g. find it on the system,
# or fetch it from the web and install it.
# See https://github.com/UCL/GreatCMakeCookOff/wiki for information

include(ChangeRootPath)

# Sets location where external project are included
if(NOT EXTERNAL_ROOT)
  set(EXTERNAL_ROOT "${CMAKE_BINARY_DIR}/external")
endif(NOT EXTERNAL_ROOT)
# Makes sure external projects are found by cmake
list(FIND CMAKE_PREFIX_PATH "${EXTERNAL_ROOT}" has_external_root)
if(has_external_root EQUAL -1)
    list(APPEND CMAKE_PREFIX_PATH "${EXTERNAL_ROOT}")
endif()

include(CachedVariables)
include(Utilities)
add_to_envvar(PKG_CONFIG_PATH "${EXTERNAL_ROOT}/lib/pkgconfig" PREPEND OS UNIX)
add_to_envvar(PKG_CONFIG_PATH
    "${EXTERNAL_ROOT}/lib64/pkgconfig" PREPEND OS UNIX)
add_to_envvar(LD_LIBRARY_PATH "${EXTERNAL_ROOT}/lib" PREPEND OS UNIX)
add_to_envvar(LD_LIBRARY_PATH "${EXTERNAL_ROOT}/lib64" PREPEND OS UNIX)
add_to_envvar(DYLD_LIBRARY_PATH "${EXTERNAL_ROOT}/lib" PREPEND OS APPLE)
add_to_envvar(DYLD_LIBRARY_PATH "${EXTERNAL_ROOT}/lib64" PREPEND OS APPLE)

# Adds a target for all external projects, so they can be made prior to others.
if(NOT TARGET lookup_dependencies)
    add_custom_target(lookup_dependencies ALL)
endif()

include(FindPackageHandleStandardArgs)
include(ExternalProject)

function(_find_lookup_recipe package OUTVAR)
    foreach(path ${CMAKE_MODULE_PATH})
      list(APPEND cmake_paths "${path}/${package}")
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

macro(_get_sane_name name OUTVAR)
    string(REGEX REPLACE "\\-" "_" ${OUTVAR} "${name}")
endmacro()

macro(_find_package_for_lookup package REQUIRED QUIET DOWNLOAD_BY_DEFAULT CHECK)
    string(TOUPPER "${package}" PACKAGE)
    set(recursive FALSE)
    _get_sane_name(${package} SANENAME)
    if(${SANENAME}_RECURSIVE) # Called by recursive step
        delete_package_variables(${package})
        set(recursive TRUE)
        if(${package}_NOFIND)
            set(${package}_FOUND TRUE CACHE BOOL "")
        endif()
    endif()
    # First try and find package (unless downloading by default)
    set(dolook TRUE)
    if(NOT recursive)
      if("${DOWNLOAD_BY_DEFAULT}" STREQUAL "TRUE")
        set(dolook FALSE)
      else()
        set(dolook TRUE)
      endif()
      if(CHECK)
          set(do_rootchange TRUE)
      endif()
    endif()
    # Figure out whether to add REQUIRED and QUIET keywords
    set(required "")
    if(recursive AND REQUIRED)
        set(required REQUIRED)
    endif()
    unset(quiet)
    if(QUIET OR (NOT QUIET AND NOT recursive))
        set(quiet QUIET)
    endif()
    if(dolook OR do_rootchange)
        if(do_rootchange)
            _set_root_path("${EXTERNAL_ROOT}" ONLY)
        else() 
            # Set external root as first place to look
            set(_SAVE_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
            list(INSERT CMAKE_PREFIX_PATH 0 "${EXTERNAL_ROOT}")
        endif()
        find_package(${package} ${ARGN} ${required} ${quiet})
        if(do_rootchange)
            _unset_root_path()
            unset(do_rootchange)
            if(${PACKAGE}_FOUND OR ${package}_FOUND)
                message(STATUS "Found ${package} in external "
                    "projects directory ${EXTERNAL_ROO}")
            endif()
        else()
            set(CMAKE_PREFIX_PATH ${_SAVE_CMAKE_PREFIX_PATH})
            if(NOT QUIET AND (${PACKAGE}_FOUND OR ${package}_FOUND))
                message(STATUS "Found ${package}")
            endif()
        endif()
    endif()
endmacro()

macro(_perform_actual_lookup package)
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
        if(NOT ${package}_QUIET)
            if(${package}_FOUND)
                message(STATUS
                    "${package} is built locally as an external project")
            elseif(NOT ${package}_DOWNLOAD_BY_DEFAULT)
                message(STATUS
                    "Will attempt to download and install ${package}")
            else()
                message(STATUS "Will download, build,"
                   " and install a local version of ${package}")
            endif()
        endif()
        set(CURRENT_LOOKUP_DIRECTORY "${${package}_LOOKUP_RECIPE_DIR}")
        if(NOT ${package}_LOOKUP_PACKAGE_INCLUSION_GUARD)
          include(${${package}_LOOKUP_RECIPE_FILE})
          set(${package}_LOOKUP_PACKAGE_INCLUSION_GUARD TRUE)
        endif()
        unset(CURRENT_LOOKUP_DIRECTORY)
        set(${package}_LOOKUP_BUILD ${${package}_KEEP} CACHE BOOL
            "Whether package is obtained from a lookup build"
            FORCE
        )
        if(TARGET ${package})
            add_dependencies(lookup_dependencies ${package})
        elseif(TARGET Lookup-${package})
            add_dependencies(lookup_dependencies Lookup-${package})
        elseif(TARGET LookUp-${package})
            add_dependencies(lookup_dependencies LookUp-${package})
        endif()
    endif()
endmacro()

# Returns name of given hook and package in hook_script variable
function(get_lookup_hookscript_name hook package)
    string(TOUPPER "${hook}" hook)
    if("${hook}" STREQUAL "POST_LOOKUP")
        set(filename "${EXTERNAL_ROOT}/hooks/post_lookup/${package}.cmake")
    elseif("${hook}" STREQUAL "INSTALL")
        set(filename "${EXTERNAL_ROOT}/hooks/install/${package}.cmake")
    else()
        set(hooks POST_LOOKUP INSTALL)
        message(FATAL_ERROR
            "Hook argument(${hook}) should be one of ${hooks}.")
    endif()
    list(LENGTH ARGN nargs)
    set(variable hook_script)
    if(NOT nargs EQUAL 0)
        list(GET ARGN 0 variable)
    endif()
    set(${variable} "${filename}" PARENT_SCOPE)
endfunction()

macro(delete_package_variables package)
    string(TOUPPER "${package}" PACKAGE)
    delete_variables(
       "LOOKUP_${package}" "LOOKUP_${PACKAGE}"
       "^${package}_LIBRAR.*" "^${PACKAGE}_LIBRAR.*"
       "^${package}_.*_LIBRAR.*" "^${PACKAGE}_.*_LIBRAR.*"
       "^${package}_INCLUD.*" "^${PACKAGE}_INCLUD.*"
       "^${package}_.*_INCLUD.*" "^${PACKAGE}_.*_INCLUD.*"
       "^${package}_FOUND" "^${PACKAGE}_FOUND"
       "^${package}_.*_FLAGS" "^${PACKAGE}_.*_FLAGS"
       "^${package}_.*CXXFLAGS" "^${PACKAGE}_.*CXXFLAGS"
       "^${package}_.*CFLAGS" "^${PACKAGE}_.*CFLAGS"
       "^${package}_.*F9.FLAGS" "^${PACKAGE}_.*F9.FLAGS"
    )
    unset(${package}_FOUND)
    unset(${PACKAGE}_FOUND)
endmacro()

# Looks for a lookup package file and includes it.
macro(lookup_package package)
    # include potential hooks
    if(${package}_BUILT_AS_EXTERNAL_PROJECT)
        get_lookup_hookscript_name(post_lookup ${package})
        if(EXISTS "${hook_script}")
            include("${hook_script}")
        endif()
    endif()

    string(TOUPPER "${package}" PACKAGE)
    cmake_parse_arguments(${package}
        "DOWNLOAD_BY_DEFAULT;REQUIRED;QUIET;KEEP;NOFIND;CHECK_EXTERNAL"
        ""
        "ARGUMENTS;COMPONENTS"
        ${ARGN}
    )
    # Set explicitly to TRUE or FALSE to simplify setting
    # ${package}_LOOKUP_BUILD
    if(${package}_KEEP)
        set(${package}_KEEP TRUE)
    else()
        set(${package}_KEEP FALSE)
    endif()
    if(NOT ${package}_REQUIRED)
      set(${package}_REQUIRED FALSE)
    endif()
    if(NOT ${package}_QUIET)
      set(${package}_QUIET FALSE)
    endif()
    if(NOT ${package}_DOWNLOAD_BY_DEFAULT)
      set(${package}_DOWNLOAD_BY_DEFAULT FALSE)
    else()
      set(${package}_DOWNLOAD_BY_DEFAULT TRUE)
    endif()
    if(NOT ${package}_CHECK_EXTERNAL)
      set(${package}_CHECK_EXTERNAL FALSE)
    endif()
    # users can request to download the package explicitly on the command-line
    if(LOOKUP_${package} OR LOOKUP_${PACKAGE})
        set(${package}_DOWNLOAD_BY_DEFAULT TRUE)
        delete_package_variables(${package})
    endif()

    # Reappends components
    if(${package}_COMPONENTS)
        list(APPEND ${package}_UNPARSED_ARGUMENTS COMPONENTS)
        list(APPEND ${package}_UNPARSED_ARGUMENTS ${${package}_COMPONENTS})
    endif()
    # First try and locate package
    _find_package_for_lookup(${package}
        ${${package}_REQUIRED} ${${package}_QUIET}
        ${${package}_DOWNLOAD_BY_DEFAULT} ${${package}_CHECK_EXTERNAL}
        ${${package}_UNPARSED_ARGUMENTS}
    )

    # Sets lower and upper case versions.
    # Otherwise some package will be registered as not found.
    # This is a problem with changing cmake practices.
    if(${PACKAGE}_FOUND AND NOT "${package}" STREQUAL "${PACKAGE}")
        set(${package}_FOUND ${${PACKAGE}_FOUND})
    endif()
    if(${package}_FOUND AND ${package}_BUILT_AS_EXTERNAL_PROJECT)
        get_lookup_hookscript_name(install ${package})
        if(EXISTS "${hook_script}")
            install(SCRIPT "${hook_script}")
        endif()
    endif()
    # If package is not found, then look for a recipe to download and build it
    if(NOT ${package}_FOUND OR ${package}_LOOKUP_BUILD OR ${package}_KEEP)
        _perform_actual_lookup(${package})
        # Sets a variable saying we are building this source externally
        set(${package}_BUILT_AS_EXTERNAL_PROJECT TRUE CACHE INTERNAL
            "${package} is the result of a local install")
    endif()
endmacro()


# Makes target depend on external dependencies
macro(depends_on_lookups TARGET)
    add_dependencies(${TARGET} lookup_dependencies)
endmacro()

# Adds an external step to an external project to rerun cmake
macro(add_recursive_cmake_step name)
    cmake_parse_arguments(recursive
        "NOCHECK" "FOUND_VAR;PACKAGE_NAME" "" ${ARGN})
    if(recursive_PACKAGE_NAME)
        set(recurse_name "${recursive_PACKAGE_NAME}")
    else()
        string(REGEX REPLACE "Look(u|U)p-?" "" recurse_name "${name}")
    endif()
    set(found_var ${recurse_name}_FOUND)
    if(recursive_FOUND_VAR)
        set(found_var ${recursive_FOUND_VAR})
    endif()

    # Only add recurse step if package not found already.
    # Once the package has been found and configured,
    # the locations and such should not change, so
    # there is no need for a recursive cmake step.
    if(NOT DEFINED ${found_var} OR NOT ${found_var})
        _get_sane_name(${recurse_name} SANENAME)
        set(cmakefile "${PROJECT_BINARY_DIR}/CMakeFiles/external")
        set(cmakefile "${cmakefile}/${name}_recursive.cmake")
        file(WRITE "${cmakefile}"
          "set(CMAKE_PROGRAM_PATH \"${EXTERNAL_ROOT}/bin\" CACHE PATH \"\")\n"
          "set(CMAKE_LIBRARY_PATH \"${EXTERNAL_ROOT}/lib\" CACHE PATH \"\")\n"
          "set(CMAKE_INCLUDE_PATH "
            "\"${EXTERNAL_ROOT}/include\" CACHE PATH \"\")\n"
          "set(CMAKE_PREFIX_PATH \"${EXTERNAL_ROOT}\" CACHE PATH \"\")\n"
          "set(${SANENAME}_RECURSIVE TRUE CACHE INTERNAL \"\")\n"
        )
        if(NOT recursive_NOCHECK)
            file(APPEND "${cmakefile}"
                "set(${SANENAME}_REQUIREDONRECURSE TRUE CACHE INTERNAL \"\")\n"
            )
        endif()
        ExternalProject_Add_Step(
            ${name} reCMake
            COMMAND
                ${CMAKE_COMMAND} -C "${cmakefile}"
                    --no-varn-unused-cli "${CMAKE_SOURCE_DIR}"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            ${recursive_UNPARSED_ARGUMENTS}
        )
        if(${SANENAME}_REQUIREDONRECURSE)
            if(NOT ${found_var} OR "${${found_var}}" STREQUAL "")
                unset(${SANENAME}_REQUIREDONRECURSE CACHE)
                message(FATAL_ERROR
                    "[${name}] Could not be downloaded and installed"
                )
            endif()
        endif()
    endif()

endmacro()

# A script that is executed once a package has been built locally
# This function should be called from the lookup recipe
function(write_lookup_hook hook package)
    cmake_parse_arguments(_wpls${package}
        "APPEND" "SCRIPTNAME;CONFIGURE" "" ${ARGN})
    if(_wpls${package}_APPEND AND _wpls${package}_CONFIGURE)
        message(FATAL_ERROR "Only one of APPEND and CONFIGURE is valid")
    endif()

    # Sets filename. That's the only difference between hooks, in practice.
    get_lookup_hookscript_name(${hook} ${package})
    # Then write or append to file.
    if(_wpls${package}_APPEND AND NOT _wpls${package}_CONFIGURE)
        file(APPEND "${hook_script}" ${_wpls${package}_UNPARSED_ARGUMENTS})
    elseif(NOT _wpls${package}_CONFIGURE)
        file(WRITE "${hook_script}" ${_wpls${package}_UNPARSED_ARGUMENTS})
    else()
        configure_file("${_wpls${package}_CONFIGURE}" "${hook_script}" @ONLY)
    endif()
    if(_wpls${package}_SCRIPTNAME)
        set(${_wpls${package}_SCRIPTNAME} "${hook_script}" PARENT_SCOPE)
    endif()
endfunction()

# Avoids anoying cmake warning, by actually using the variables.
# The will be if the appropriate find_* is used. But won't be otherwise.
if(CMAKE_PROGRAM_PATH)
endif()
if(CMAKE_LIBRARY_PATH)
endif()
if(CMAKE_INCLUDE_PATH)
endif()
