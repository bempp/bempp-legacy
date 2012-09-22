configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
	${CMAKE_BINARY_DIR}/include/dune/foamgrid/config.h)

include_directories(${CMAKE_BINARY_DIR}/include/dune/foamgrid)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

