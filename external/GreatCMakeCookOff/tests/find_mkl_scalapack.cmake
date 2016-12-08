cmake_minimum_required(VERSION 2.8)
project(mkl_scalapack C)

find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()

find_package(MKL REQUIRED COMPONENTS ScaLAPACK)
find_package(Threads REQUIRED)
find_package(MPI REQUIRED)

add_definitions(${MKL_DEFINITIONS})
include_directories(${MKL_INCLUDE_DIRS})

message(STATUS "MKL_DEFINITIONS = ${MKL_DEFINITIONS}")
message(STATUS "MKL_INCLUDE_DIRS = ${MKL_INCLUDE_DIRS}")
message(STATUS "MKL_LIBRARIES = ${MKL_LIBRARIES}")

file(WRITE "${CMAKE_SOURCE_DIR}/mkl_scalapack.c"
    "#include <stddef.h>\n"
    "#include <mkl.h>\n"
    "extern void blacs_setup_(int *, int *);"
    "extern void pdgemr2d_(int *, int *, double *, int *, int *, int *, double *,\n"
    "int *, int *, int *, int *);\n"
    "int main() {\n"
    "  blacs_setup_(NULL, NULL);\n"
    "  pdgemr2d_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);\n"
    "  return 0;\n"
    "}\n"
)

add_executable(mkl_scalapack mkl_scalapack.c)
target_link_libraries(mkl_scalapack ${MKL_LIBRARIES} ${MPI_C_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Using GCC requires -lm
include(CheckFunctionExists)

CHECK_FUNCTION_EXISTS(logf LOGF_EXISTS)
if(NOT LOGF_EXISTS)
  unset(LOGF_EXISTS)
  unset(LOGF_EXISTS CACHE)
  list(APPEND CMAKE_REQUIRED_LIBRARIES m)
  CHECK_FUNCTION_EXISTS(logf LOGF_EXISTS)
  if(LOGF_EXISTS)
    target_link_libraries(mkl_scalapack m)
  else()
    message(FATAL_ERROR "No logf() found")
  endif()
endif()
