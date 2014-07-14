if (SWIG_FOUND)
	file(GLOB_RECURSE PYTHON_EXAMPLE_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
	if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    	foreach (f ${PYTHON_EXAMPLE_SOURCES})
        	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${f}
                       	${CMAKE_CURRENT_BINARY_DIR}/${f}
                       	COPYONLY)
    	endforeach ()
	endif ()
install(FILES ${PYTHON_EXAMPLE_SOURCES} DESTINATION ${SHARE_INSTALL_PATH}/bempp/examples/python)
endif ()

