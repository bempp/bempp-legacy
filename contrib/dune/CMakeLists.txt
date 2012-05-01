set(ALUGRID_ROOT_DIR "" CACHE PATH "Full path to the ALUGrid root directory")
if (WITH_ALUGRID)
        # Variable used in dunecontrol.conf.in
        set(DUNE_WITH_ALUGRID_SETTING "--with-alugrid=${ALUGRID_ROOT_DIR}" CACHE INTERNAL "")
else ()
        set(DUNE_WITH_ALUGRID_SETTING "" CACHE INTERNAL "")
endif ()
configure_file(dunecontrol.conf.in dunecontrol.conf)

add_custom_target(dune
	DEPENDS dune-grid dune-common dune-localfunctions
)

include(ExternalProject)
ExternalProject_Add(
	dune-grid
	PREFIX ${CMAKE_SOURCE_DIR}/dune
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/dune/dune-grid
	BINARY_DIR ${CMAKE_SOURCE_DIR}/dune/dune-grid
	DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/dune
	URL http://www.dune-project.org/download/2.1/dune-grid-2.1.0.tar.gz
	CONFIGURE_COMMAND ""
	PATCH_COMMAND patch ${CMAKE_SOURCE_DIR}/dune/dune-grid/dune/grid/utility/structuredgridfactory.hh < ${CMAKE_SOURCE_DIR}/patch/dune-grid.patch
	BUILD_COMMAND ""
	UPDATE_COMMAND ""
	INSTALL_COMMAND ""
)

ExternalProject_Add(
	dune-localfunctions
	PREFIX ${CMAKE_SOURCE_DIR}/dune
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/dune/dune-localfunctions
	BINARY_DIR ${CMAKE_SOURCE_DIR}/dune/dune-localfunctions
	DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/dune
	URL http://www.dune-project.org/download/2.0/dune-localfunctions-2.0.tar.gz
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
	UPDATE_COMMAND ""
	INSTALL_COMMAND ""
)

ExternalProject_Add(
	dune-common
	DEPENDS dune-grid dune-localfunctions
	PREFIX ${CMAKE_SOURCE_DIR}/dune
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/dune/dune-common
	BINARY_DIR ${CMAKE_SOURCE_DIR}/dune/dune-common
	DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/dune
	URL http://www.dune-project.org/download/2.1/dune-common-2.1.0.tar.gz
	CONFIGURE_COMMAND ""
#	BUILD_COMMAND DUNE_CONTROL_PATH=${CMAKE_SOURCE_DIR}/dune CONFIGURE_FLAGS='--enable-shared --disable-static' ${CMAKE_SOURCE_DIR}/dune/dune-common/bin/dunecontrol all
	BUILD_COMMAND DUNE_CONTROL_PATH=${CMAKE_SOURCE_DIR}/dune ${CMAKE_SOURCE_DIR}/dune/dune-common/bin/dunecontrol --opts=${CMAKE_SOURCE_DIR}/dune/dunecontrol.conf all 
	UPDATE_COMMAND ""
	INSTALL_COMMAND ""
)	

