cmake_minimum_required (VERSION 2.8)
project (contrib)

# Set Cache Variables

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
    "${CMAKE_SOURCE_DIR}" CACHE PATH "Install prefix" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)



# Add directories 

if (WITH_DUNE)
   add_subdirectory(dune)
endif (WITH_DUNE)

