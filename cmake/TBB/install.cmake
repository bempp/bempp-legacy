# simple CMake file to install tbb from its extracted download file to places we want.
cmake_minimum_required(VERSION 2.8)

# Find patch files and apply
file(GLOB PATCH_FILES ${ROOT}/cmake/patches/tbb/*.patch)
find_program(PATCH_EXECUTABLE patch)
foreach(patch ${PATCH_FILES})
  execute_process(
      COMMAND ${PATCH_EXECUTABLE} -N < ${patch}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
endforeach()

include(DetectArchitecture)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(libglob "${PROJECT_SOURCE_DIR}/lib/lib*.*")
elseif("${DETECTED_ARCHITECTURE}" STREQUAL "x86_64")
    set(libglob "${PROJECT_SOURCE_DIR}/lib/intel64/gcc4.4/lib*.*")
elseif("${DETECTED_ARCHITECTURE}" STREQUAL "i386")
    set(libglob "${PROJECT_SOURCE_DIR}/lib/ia32/gcc4.4/lib*.*")
elseif("${DETECTED_ARCHITECTURE}" STREQUAL "ia64")
    set(libglob "${PROJECT_SOURCE_DIR}/lib/ia64/gcc4.4/lib*.*")
else()
    message(FATAL_ERROR "Could not detect architecture")
endif()
# install libraries
file(GLOB LIBRARIES "${libglob}")
install(FILES ${LIBRARIES} DESTINATION lib)

# Now deal with recursive install of headers
function(add_headers directory relative)
    # Find headers and install them
    file(GLOB HEADERS RELATIVE "${relative}" "${directory}/*.h")
    if(HEADERS)
        list(GET HEADERS 0 header_file)
        get_filename_component(relative_directory "${header_file}" PATH)
        install(FILES ${HEADERS} DESTINATION "${relative_directory}")
    endif()
    # Recurse
    file(GLOB subdirectories "${directory}/*/")
    foreach(subdirectory ${subdirectories})
        if(IS_DIRECTORY "${subdirectory}")
            add_headers("${subdirectory}" "${relative}")
        endif()
    endforeach()
endfunction()

add_headers("${PROJECT_SOURCE_DIR}/include/tbb" "${PROJECT_SOURCE_DIR}")
