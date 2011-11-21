# Please observe: line 32 and 33 only work with newer DUNE modules.
# IF you use an old module, that has a symlink "dune -> .", you have to remove "/dune/" from these lines!

MACRO( ADD_CXX_FLAGS )
  FOREACH( ARG ${ARGN} )
    ADD_DEFINITIONS( ${ARG} )
    LIST( APPEND MY_CXX_FLAGS ${ARG} )
  ENDFOREACH( ARG )
ENDMACRO( ADD_CXX_FLAGS )

MACRO( INCLUDE_DIR )
  FOREACH( ARG ${ARGN} )
    INCLUDE_DIRECTORIES( ${ARG} )
    LIST( APPEND MY_CXX_FLAGS "-I${ARGN}" )
    ENDFOREACH( ARG )
ENDMACRO( INCLUDE_DIR )

MACRO( HEADERCHECK )
  ADD_CUSTOM_TARGET( headercheck )
  FOREACH( HEADER ${ARGN} )
    GET_FILENAME_COMPONENT( fn ${HEADER} NAME )
    SET( TEST_NAME "headercheck_${fn}" )
    ADD_CUSTOM_TARGET( ${TEST_NAME} ${CMAKE_CXX_COMPILER} ${MY_CXX_FLAGS} -c ${HEADER} -o ${CMAKE_BINARY_DIR}/tests/${TEST_NAME}.o )
    ADD_TEST( ${TEST_NAME} ${CMAKE_CXX_COMPILER} ${MY_CXX_FLAGS} -c ${HEADER} -o ${CMAKE_BINARY_DIR}/tests/${TEST_NAME}.o )
    add_dependencies( headercheck ${TEST_NAME} )
  ENDFOREACH( HEADER )
ENDMACRO( HEADERCHECK )

MACRO( ADD_DUNE_MODULES )
  FOREACH( MODULE ${ARGN} )
	  INCLUDE_DIR( ${DUNE_PATH}/dune-${MODULE} )
	  LINK_DIRECTORIES( ${DUNE_PATH}/dune-${MODULE}/dune/${MODULE}/.libs )
	  FILE( GLOB_RECURSE tmp_header "${DUNE_PATH}/dune-${MODULE}/dune/${MODULE}/*.hh" )
    LIST( APPEND DUNE_HEADERS ${tmp_header} )
  ENDFOREACH( MODULE )
ENDMACRO( ADD_DUNE_MODULES )

