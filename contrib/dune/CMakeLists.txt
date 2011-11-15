add_custom_target(dune
	DEPENDS dune-grid dune-foamgrid dune-common
)

include(ExternalProject)
ExternalProject_Add(
	dune-grid
	PREFIX ${CMAKE_BINARY_DIR}/dune
	SOURCE_DIR ${CMAKE_BINARY_DIR}/dune/dune-grid
	BINARY_DIR ${CMAKE_BINARY_DIR}/dune/dune-grid
	DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/dune
	URL http://www.dune-project.org/download/2.1/dune-grid-2.1.0.tar.gz
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
	UPDATE_COMMAND ""
	INSTALL_COMMAND ""
)

ExternalProject_Add(
	dune-foamgrid
	PREFIX ${CMAKE_BINARY_DIR}/dune
	SOURCE_DIR ${CMAKE_BINARY_DIR}/dune/dune-foamgrid
	BINARY_DIR ${CMAKE_BINARY_DIR}/dune/dune-foamgrid
	DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/dune
	CONFIGURE_COMMAND ""
	URL file://${CMAKE_SOURCE_DIR}/contrib/dune/dune-foamgrid-2.1.tar.gz
	BUILD_COMMAND ""
	UPDATE_COMMAND ""
	INSTALL_COMMAND ""
)



ExternalProject_Add(
	dune-common
	DEPENDS dune-grid dune-foamgrid
	PREFIX ${CMAKE_BINARY_DIR}/dune
	SOURCE_DIR ${CMAKE_BINARY_DIR}/dune/dune-common
	BINARY_DIR ${CMAKE_BINARY_DIR}/dune/dune-common
	DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/dune
	URL http://www.dune-project.org/download/2.1/dune-common-2.1.0.tar.gz
	CONFIGURE_COMMAND ""
	BUILD_COMMAND DUNE_CONTROL_PATH=${CMAKE_BINARY_DIR}/dune ${CMAKE_BINARY_DIR}/dune/dune-common/bin/dunecontrol all
	UPDATE_COMMAND ""
	INSTALL_COMMAND ""
)	

