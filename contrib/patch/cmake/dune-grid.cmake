configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
	${CMAKE_BINARY_DIR}/include/dune/grid/config.h)

include_directories(${CMAKE_BINARY_DIR}/include/dune/grid)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB QUADRATURERULES_LIB_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/dune/grid/common/quadraturerules/*.cc)
file(GLOB ONEDGRID_LIB_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/dune/grid/onedgrid/*.cc)

set(DGFPARSER_LIB_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/dune/grid/io/file/dgfparser/dgfparserblocks.cc
	${CMAKE_CURRENT_SOURCE_DIR}/dune/grid/io/file/dgfparser/dgfprojectionblock.cc
	${CMAKE_CURRENT_SOURCE_DIR}/dune/grid/io/file/dgfparser/dgfparser.cc
)	

add_library(dunegrid SHARED 
	${QUADRATURERULES_LIB_SOURCES}
	${ONEDGRID_LIB_SOURCES}
	${DGFPARSER_LIB_SOURCES}
		)


#set_target_properties(quadraturerules PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(dunegrid dunecommon)

install(TARGETS dunegrid 
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

