set(
  OPENCL_INC_SEARCH_PATH
  ${OPENCL_INCLUDE_DIR}
  ${OPENCL_DIR}/include
  ${OPENCL_DIR}/OpenCL/common/inc
  $ENV{OPENCL_INCLUDE_DIR}
  $ENV{OPENCL_DIR}/include
  $ENV{OPENCL_DIR}/OpenCL/common/inc
  /usr/local/include
  /usr/include
  /usr/local/cuda/include
  # Append additional search directories for OpenCL include files here.
)

set(
  OPENCL_LIB_SEARCH_PATH
  ${OPENCL_LIBRARY_DIR}
  ${OPENCL_DIR}/lib
  ${OPENCL_DIR}/lib/x86
  $ENV{OPENCL_LIBRARY_DIR}
  $ENV{OPENCL_DIR}/lib
  $ENV{OPENCL_DIR}/lib/x86
  /usr/local/lib64
  /usr/local/lib
  /usr/lib64
  /usr/lib
  # Append additional search directories for OpenCL libraries here.
)

if("${CMAKE_SYSTEM_PROCESSOR}" EQUAL "x86_64")
  set(
        OPENCL_LIB_SEARCH_PATH
        ${OPENCL_LIB_SEARCH_PATH}
        ${OPENCL_DIR}/OpenCL/common/lib/Linux64
        $ENV{OPENCL_DIR}/OpenCL/common/lib/Linux64
        )
else("${CMAKE_SYSTEM_PROCESSOR}" EQUAL "x86_64")
  set(
        OPENCL_LIB_SEARCH_PATH
        ${OPENCL_LIB_SEARCH_PATH}
        ${OPENCL_DIR}/OpenCL/common/lib/Linux32
        $ENV{OPENCL_DIR}/OpenCL/common/lib/Linux32
        )
endif("${CMAKE_SYSTEM_PROCESSOR}" EQUAL "x86_64")

# Do not change anything of below (except if bugs encountered).
find_path(
  OPENCL_INCLUDE_DIR
  NAMES CL/cl.h
  PATHS ${OPENCL_INC_SEARCH_PATH}
  )

find_library(
  OPENCL_LIBRARY
  NAMES OpenCL
  PATHS ${OPENCL_LIB_SEARCH_PATH}
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  OPENCL
  DEFAULT_MSG
  OPENCL_LIBRARY OPENCL_INCLUDE_DIR
  )

if(OPENCL_FOUND)
  set(OPENCL_LIBRARIES ${OPENCL_LIBRARY})
else(OPENCL_FOUND)
  set(OPENCL_LIBRARIES)
endif(OPENCL_FOUND)

mark_as_advanced(
  OPENCL_INCLUDE_DIR
  OPENCL_LIBRARY
  )