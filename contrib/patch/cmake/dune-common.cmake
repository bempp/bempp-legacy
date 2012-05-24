configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
	${CMAKE_BINARY_DIR}/include/dune/common/config.h)

include_directories(${CMAKE_BINARY_DIR}/include/dune/common)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(dunecommon SHARED 
	${CMAKE_CURRENT_SOURCE_DIR}/dune/common/configparser.cc
	${CMAKE_CURRENT_SOURCE_DIR}/dune/common/exceptions.cc
	${CMAKE_CURRENT_SOURCE_DIR}/dune/common/ios_state.cc
	${CMAKE_CURRENT_SOURCE_DIR}/dune/common/parametertree.cc
	${CMAKE_CURRENT_SOURCE_DIR}/dune/common/parametertreeparser.cc
	${CMAKE_CURRENT_SOURCE_DIR}/dune/common/path.cc
	${CMAKE_CURRENT_SOURCE_DIR}/dune/common/stdstreams.cc
	)
set_target_properties(dunecommon PROPERTIES LINKER_LANGUAGE CXX)

install(TARGETS dunecommon
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

