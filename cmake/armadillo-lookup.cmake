# Looks for armadillo. If not found, download and install it.
find_package(Armadillo)
if(NOT ARMADILLO_INCLUDE_DIR)
  message(STATUS "[Armadillo] not found. Will attempt to download and install it.")
  include(ExternalProject)
  
  # Downloads armadillo
  ExternalProject_Add(
      Armadillo
      PREFIX ${EXTERNAL_ROOT}
      URL http://sourceforge.net/projects/arma/files/armadillo-4.000.4.tar.gz
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT}
      TIMEOUT 10
      # Wrap download, configure and build steps in a script to log output
      LOG_DOWNLOAD ON
      LOG_CONFIGURE ON
      LOG_BUILD ON
  )
  # Rerun cmake to capture new armadillo install
  add_recursive_cmake_step(Armadillo Armadillo_FOUND DEPENDEES install)
  add_compile_options(ARMA_USE_LAPACK)
  add_compile_options(ARMA_USE_BLAS)
endif()
