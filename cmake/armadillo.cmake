# Look for boost the easy way first
if(NOT ARMADILLO_INCLUDE_DIR)

  # Not clear that parsing works in find_package. So setting variable directly.
  find_package(Armadillo)
  if(Armadillo_FOUND)
    set(
      ARMADILLO_INCLUDE_DIR ${Armadillo_INCLUDE_DIR}
      CACHE INTERNAL
      "Path to armadillo include directory"
    )
  endif()
endif()


# Install it if not found
if(NOT Armadillo_FOUND)
  if(NOT EXTERNAL_ROOT)
    set(EXTERNAL_ROOT ${CMAKE_BINARY_DIR}/external)
  endif(NOT EXTERNAL_ROOT)

  message(STATUS "Armadillo not found. Will attempt to download it.")
  include(ExternalProject)
  
  # Downloads armadillo
  ExternalProject_Add(
      armadillo
      PREFIX ${EXTERNAL_ROOT}
      URL http://sourceforge.net/projects/arma/files/armadillo-4.000.4.tar.gz
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT}
      TIMEOUT 10
      # Wrap download, configure and build steps in a script to log output
      LOG_DOWNLOAD ON
      LOG_CONFIGURE ON
      LOG_BUILD ON
  )
  set(
    ARMADILLO_INCLUDE_DIR ${EXTERNAL_ROOT}/include
    CACHE INTERNAL
    "Path to armadillo include directory"
  )
endif()

