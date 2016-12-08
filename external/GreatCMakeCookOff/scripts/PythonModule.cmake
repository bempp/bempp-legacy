include(FindPackageHandleStandardArgs)
include(PythonInstall)
include(TargetCopyFiles)
include(EnvironmentScript)
include(FilterList)

if(NOT PYTHON_BINARY_DIR)
    set(PYTHON_BINARY_DIR "${PROJECT_BINARY_DIR}/python_package"
        CACHE PATH "Location of python package in build tree")
endif()
add_to_python_path(${PYTHON_BINARY_DIR})
set(DEPS_SCRIPT
    ${CMAKE_CURRENT_LIST_DIR}/find_cython_deps.py
    CACHE INTERNAL "Script to determine cython dependencies"
)

function(_pm_location_and_name module module_LOCATION)
    string(REGEX REPLACE "\\." "/" location "${module}")
    get_filename_component(submodule "${location}" NAME)

    set(submodule ${submodule} PARENT_SCOPE)

    if(NOT "${module_LOCATION}" STREQUAL "")
        set(location ${module_LOCATION} PARENT_SCOPE)
    endif()

    set(location ${location} PARENT_SCOPE)
endfunction()

function(_pm_default)
    if(${module}_NOINSTALL)
        set(do_install FALSE PARENT_SCOPE)
    else()
        set(do_install TRUE PARENT_SCOPE)
    endif()
    if(NOT ${module}_CPP)
        set(${module}_CPP "" PARENT_SCOPE)
    else()
        set(${module}_CPP CPP PARENT_SCOPE)
    endif()

    unset(excluded)
    if(${module}_EXCLUDE)
        file(GLOB files RELATIVE
            "${CMAKE_CURRENT_SOURCE_DIR}" ${${module}_EXCLUDE})
        list(APPEND excluded ${files})
    endif()

    unset(sources)
    if(NOT "${${module}_GLOB}" STREQUAL "")
        file(GLOB sources
            RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" ${${module}_GLOB})
    endif()
    list(APPEND sources ${${module}_SOURCES} ${${module}_UNPARSED_ARGUMENTS})
    list(REMOVE_DUPLICATES sources)
    if(NOT "${excluded}" STREQUAL "")
        list(REMOVE_ITEM sources ${excluded})
        file(GLOB patterns ${excluded})
        if(NOT "${patterns}" STREQUAL "")
            list(REMOVE_ITEM sources ${patterns})
        endif()
    endif()

    if("${sources}" STREQUAL "")
        message(FATAL_ERROR "Python module has no sources")
    endif()
    set(ALL_SOURCES ${sources} PARENT_SCOPE)
endfunction()

function(get_pyx_dependencies SOURCE OUTVAR)
    set(local_python "${LOCAL_PYTHON_EXECUTABLE}")
    if(NOT "${local_python}")
        set(local_python ${PYTHON_EXECUTABLE})
    endif()
    execute_process(
        COMMAND ${local_python} ${DEPS_SCRIPT} ${SOURCE} ${ARGN}
        RESULT_VARIABLE RESULT
        OUTPUT_VARIABLE OUTPUT
        ERROR_VARIABLE ERROR
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    if("${RESULT}" STREQUAL "0")
        set(${OUTVAR} ${OUTPUT} PARENT_SCOPE)
    else()
        message("Error: ${ERROR}")
        message("Output: ${OUTPUT}")
        message(FATAL_ERROR "Error while computing cython dependencies")
    endif()
endfunction()

function(_pm_add_fake_init location)
    set(fake_init_file "${PYTHON_BINARY_DIR}/${location}/__init__.py")
    if(NOT EXISTS "${fake_init_file}")
        file(WRITE "${fake_init_file}" "# Empty file added by CMake")
    endif()
    if(do_install)
        install_python(FILES "${fake_init_file}" DESTINATION ${location})
    endif()
endfunction()

function(python_extension_targetname outvar module)
    cmake_parse_arguments(_pm_tgname
        ""
        "MODULE_TARGET"
        ";"
        ${ARGN}
    )
    set(module_target ${module}-ext)
    if(NOT "${_pm_tgname_MODULE_TARGET}" STREQUAL "")
        set(module_target ${_pm_tgname_MODULE_TARGET})
    endif()
    set(${outvar} ${module_target} PARENT_SCOPE)
endfunction()

function(_pm_add_python_extension module)
    string(REGEX REPLACE "/" "_" ext "ext.${module}")
    cmake_parse_arguments(${ext}
        ""
        "INSTALL;TARGET;LOCATION;EXTENSION;MODULE_TARGET"
        "SOURCES;LIBRARIES;DEPENDENCIES;OBJECTS"
        ${ARGN}
    )
    if("${${ext}_SOURCES}" STREQUAL "")
        return()
    endif()

    include_directories(${PYTHON_INCLUDE_DIRS})
    if(NUMPY_INCLUDE_DIRS)
        include_directories(${NUMPY_INCLUDE_DIRS})
    endif()

    set(location ${${ext}_LOCATION})
    set(container_target ${${ext}_TARGET})
    python_extension_targetname(module_target
        ${${ext}_TARGET} MODULE_TARGET ${${ext}_MODULE_TARGET})

    add_library(${module_target} MODULE ${${ext}_SOURCES} ${${ext}_OBJECTS})
    target_link_libraries(${module_target} ${PYTHON_LIBRARIES})
    set(output_dir "${location}")
    if(NOT IS_ABSOLUTE "${location}")
        set(output_dir "${PYTHON_BINARY_DIR}/${location}")
    endif()
    set_target_properties(${module_target}
        PROPERTIES
        OUTPUT_NAME "${${ext}_EXTENSION}"
        PREFIX "" SUFFIX ".so"
        LIBRARY_OUTPUT_DIRECTORY "${output_dir}"
    )
    if(${ext}_LIBRARIES)
        target_link_libraries(${module_target} ${${ext}_LIBRARIES})
    endif()
    add_dependencies(${container_target} ${module_target})
    if(NOT "${${ext}_DEPENDENCIES}" STREQUAL "")
        add_dependencies(${module_target} ${${ext}_DEPENDENCIES})
    endif()

    if(${${ext}_INSTALL})
        install_python(TARGETS ${module_target} DESTINATION "${location}")
    endif()
endfunction()

function(_pm_add_pure_python)
    string(REGEX REPLACE "/" "_" py "py.${module}")
    cmake_parse_arguments(${py}
        ""
        "INSTALL;TARGET;LOCATION"
        "SOURCES"
        ${ARGN}
    )
    if("${${py}_SOURCES}" STREQUAL "")
        return()
    endif()

    file(RELATIVE_PATH targetname_copy "${PROJECT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}")
    string(REGEX REPLACE "( |/)" "_" targetname_copy "${targetname_copy}")
    if(targetname_copy STREQUAL "")
        set(targetname_copy "${${py}_TARGET}-copy")
    else()
        set(targetname_copy "${targetname_copy}-copy")
    endif()

    add_copy_files(${${py}_TARGET}
        FILES ${${py}_SOURCES}
        DESTINATION "${PYTHON_BINARY_DIR}/${${py}_LOCATION}"
    )
    if(${${py}_INSTALL})
        install_python(FILES ${${py}_SOURCES} DESTINATION ${${py}_LOCATION})
    endif()
endfunction()

function(_pm_add_headers module)
    string(REGEX REPLACE "/" "_" h "h.${module}")
    cmake_parse_arguments(${h}
        ""
        "INSTALL;LOCATION;DESTINATION;TARGET"
        "SOURCES"
        ${ARGN}
    )
    if(NOT ${${h}_INSTALL})
        return()
    endif()

    set(headers ${${h}_SOURCES})
    if("${headers}" STREQUAL "")
        return()
    endif()
    list(REMOVE_DUPLICATES headers)

    string(FIND "${module}" "." first_dot)
    if(first_dot EQUAL -1)
        set(base_module ${module})
    else()
        string(SUBSTRING "${module}" 0 ${first_dot} base_module)
    endif()

    set(header_destination ${base_module}/include/${${h}_LOCATION})
    if(${h}_DESTINATION)
        string(REGEX REPLACE "\\." "/" header_destination ${${h}_DESTINATION})
    endif()

    if(NOT IS_ABSOLUTE "${header_destination}")
        add_copy_files(${${h}_TARGET}
            FILES ${headers}
            DESTINATION "${PYTHON_BINARY_DIR}/${header_destination}"
        )
    endif()
    install_python(FILES ${headers}
        DESTINATION ${header_destination}
        COMPONENT dev
    )
endfunction()

function(_pm_add_cythons module)
    string(REGEX REPLACE "/" "_" cys "cys.${module}")
    cmake_parse_arguments(${cys} "" "" "SOURCES" ${ARGN})
    foreach(source ${${cys}_SOURCES})
        _pm_add_cython(${module} ${source} ${${cys}_UNPARSED_ARGUMENTS})
    endforeach()
endfunction()

function(_pm_add_cython module source)

    string(REGEX REPLACE "/" "_" cy "cy.${module}")
    cmake_parse_arguments(${cy} "CPP" "TARGET;STARTLINE" "" ${ARGN})
    # Creates command-line arguments for cython for include directories
    get_property(included_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
    set(inclusion)
    foreach(included ${included_dirs})
        set(inclusion ${inclusion} -I${included})
    endforeach()

    # Computes dependencies
    get_pyx_dependencies(${source} DEPENDENCIES ${included_dirs})

    # Call cython
    string(REGEX REPLACE "\\.pyx" "" cy_module "${source}")
    string(REGEX REPLACE "\\.pxd" "" cy_module "${cy_module}")
    string(REGEX REPLACE "/" "." cy_module "${cy_module}")
    unset(arguments)
    if(cython_EXECUTABLE)
        set(arguments ${cython_EXECUTABLE})
    elseif(LOCAL_PYTHON_EXECUTABLE)
        set(arguments ${LOCAL_PYTHON_EXECUTABLE} -m cython)
    else()
        set(arguments ${PYTHON_EXECUTABLE} -m cython)
    endif()
    if(${cy}_CPP)
        set(generated_source "cython_${cy_module}.cc")
        list(APPEND arguments --cplus)
    else()
        set(generated_source "cython_${cy_module}.c")
    endif()

    # Create C source from cython
    list(APPEND arguments
        "${source}"
        -o "${CMAKE_CURRENT_BINARY_DIR}/${generated_source}" ${inclusion}
    )
    if(NOT CMAKE_VERSION VERSION_LESS "2.8.10")
        set(cond "$<$<OR:$<CONFIG:RelWithDebInfo>,$<CONFIG:Debug>>:")
        set(cmdline
            "${cond}--gdb> ${cond}--gdb-outdir> ${cond}${PROJECT_BINARY_DIR}/cython_debug_files>")
    endif()
    add_custom_command(
        OUTPUT "${generated_source}"
        COMMAND ${arguments}
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        DEPENDS ${DEPENDENCIES}
        COMMENT "Generating c/c++ source ${source} with cython (${directory})"
    )
    # Adds line at the top of the generated file
    # Its necessary to deal with some ctype.h vs cctype unpleasantness
    if(${cy}_STARTLINE)
      get_filename_component(directory "${generated_source}" DIRECTORY)
      get_filename_component(filename "${generated_source}" NAME)
      if("${directory}" STREQUAL "")
          set(c_source "startlines.${filename}")
      else()
          set(c_source "${ddir}/startlines.${filename}")
      endif()
      if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/prepend_line.sh")
          file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/prepend_line.sh"
            "cd ${CMAKE_CURRENT_BINARY_DIR};"
            "echo \"${${cy}_STARTLINE}\" > $2;"
            "cat $1 >> $2;"
          )
      endif()
      add_custom_command(
          OUTPUT "${c_source}"
          COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/prepend_line.sh ${generated_source} ${c_source}
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
          DEPENDS "${generated_source}"
          COMMENT "Adding starting line to generated cython source ${source} (${directory})"
      )
    else()
      set(c_source "${generated_source}")
    endif()

    # Extension name
    string(REGEX MATCH "[^\\.]*$" extension ${cy_module})
    if("${extension}" STREQUAL "")
        set(extension ${cy_module})
    endif()

    _pm_cython_full_module_name(full_module ${module} "${source}")
    cython_extension_targetname(targetname ${module} "${source}")
    # Add python module
    _pm_add_python_extension(${full_module}
        TARGET ${${cy}_TARGET}
        MODULE_TARGET ${targetname}
        EXTENSION ${extension}
        SOURCES ${c_source}
        ${${cy}_UNPARSED_ARGUMENTS}
    )
endfunction()

function(_pm_cython_full_module_name outvar module source)
    string(REGEX REPLACE "\\.pyx" "" cy_module "${source}")
    string(REGEX REPLACE "\\.pxd" "" cy_module "${cy_module}")
    string(REGEX REPLACE "\\.py" "" cy_module "${cy_module}")
    get_filename_component(cy_module "${cy_module}" NAME)
    if("${cy_module}" MATCHES "^${module}")
        set(full_module ${cy_module})
    else()
        set(full_module ${module}.${cy_module})
    endif()
    set(${outvar} ${full_module} PARENT_SCOPE)
endfunction()

function(cython_extension_targetname outvar module source)
    _pm_cython_full_module_name(full_module ${module} "${source}")
    python_extension_targetname(targetname ${module}
        MODULE_TARGET ${full_module}-cython)
    set(${outvar} ${targetname} PARENT_SCOPE)
endfunction()

function(_pm_get_confed_filename filename OUTPUT)
    get_filename_component(filename "${filename}" ABSOLUTE)
    file(RELATIVE_PATH relfile "${CMAKE_CURRENT_SOURCE_DIR}" "${filename}")
    if("${relfile}" MATCHES "\\.\\./")
        file(RELATIVE_PATH relfile "${CMAKE_CURRENT_BINARY_DIR}" "${filename}")
        if("${relfile}" MATCHES "\\.\\./")
            message(FATAL_ERROR "File ${filename} is not in build or source "
                "directory or subdirectory.")
        endif()
    endif()
    set(${OUTPUT} "${CMAKE_CURRENT_BINARY_DIR}/${relfile}" PARENT_SCOPE)
endfunction()

function(add_python_module module)

    # Parses arguments
    cmake_parse_arguments(${module}
        "FAKE_INIT;NOINSTALL;INSTALL;CPP"
        "HEADER_DESTINATION;TARGETNAME;LOCATION;OUTPUT_PYTHON_SOURCES"
        "SOURCES;EXCLUDE;LIBRARIES;GLOB;OBJECTS"
        ${ARGN}
    )
    # Sets submodule, location, and module from module
    _pm_location_and_name(${module} "${${module}_LOCATION}")

    # Sets defaults, do_install, and  ALL_SOURCES
    _pm_default()
    set(targetname ${module})
    if(${module}_TARGETNAME)
        set(targetname ${${module}_TARGETNAME})
    endif()
    # creates a global target
    if(NOT TARGET ${targetname})
        add_custom_target(${targetname} ALL)
    endif()

    # Figures out C/C++/HEADERS/Python sources
    filter_list(C_SOURCES ALL_SOURCES ".*\\.c$")
    filter_list(C_HEADERS ALL_SOURCES ".*\\.h$")
    filter_list(CPP_SOURCES ALL_SOURCES ".*\\.cpp$" ".*\\.cc$")
    filter_list(CPP_HEADERS ALL_SOURCES ".*\\.hpp" ".*\\.h")
    filter_list(PY_SOURCES ALL_SOURCES ".*\\.py$")
    filter_list(CY_SOURCES ALL_SOURCES ".*\\.pyx")
    filter_list(CY_HEADERS ALL_SOURCES ".*\\.pxd")

    if(C_SOURCES OR CPP_SOURCES)
        if(PY_SOURCES OR CY_SOURCES)
            message(FATAL_ERROR "Python/Cython and C sources in same call"
                " to add_python_module.\n"
                "Please split into separate pure C extensions from othes."
            )
        endif()
    endif()

    # Now for the actual meat

    # First adds fake init if necessary
    if(${module}_FAKE_INIT)
        if(C_SOURCES OR CPP_SOURCES)
            message(FATAL_ERROR
                "FAKE_INIT AND C/C++ extensions are incompatible")
        endif()
        _pm_add_fake_init(${location})
    endif()


    # Then compiles an extension if C/C++ sources
    get_filename_component(extension_location "${location}" PATH)
    _pm_add_python_extension(${module}
        TARGET ${targetname}
        INSTALL ${do_install}
        EXTENSION ${submodule}
        LOCATION ${extension_location}
        LIBRARIES ${${module}_LIBRARIES}
        SOURCES ${C_SOURCES} ${CPP_SOURCES}
        OBJECTS ${${module}_OBJECTS}
    )

    # Then copy/install pure python files
    _pm_add_pure_python(${module}
        TARGET ${targetname}
        INSTALL ${do_install}
        LOCATION ${location}
        SOURCES ${PY_SOURCES}
    )

    # Then copy/install header files
    _pm_add_headers(${module}
        TARGET ${targetname}
        LOCATION ${location}
        DESTINATION ${${module}_HEADER_DESTINATION}
        SOURCES ${CPP_HEADERS} ${C_HEADERS} ${CY_HEADERS}
        INSTALL ${do_install}
    )

    # Then create cython extensions
    _pm_add_cythons(${module}
        ${${module}_CPP}
        LOCATION ${location}
        INSTALL ${do_install}
        LIBRARIES ${${module}_LIBRARIES}
        TARGET ${targetname}
        OBJECTS ${${module}_OBJECTS}
        SOURCES ${CY_SOURCES}
    )

    # Outputs pure python sources if requested.
    # This is used mainly by add_pytest. It makes configurable tests trivial to
    # add.
    if(NOT ${${module}_OUTPUT_PYTHON_SOURCES} STREQUAL "")
        set(${${module}_OUTPUT_PYTHON_SOURCES} ${PY_SOURCES} PARENT_SCOPE)
    endif()

endfunction()
