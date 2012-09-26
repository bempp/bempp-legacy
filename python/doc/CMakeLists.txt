#Configure All Option files

# Add a python_doc target to generate Python API documentation with Sphinx
find_package(Sphinx)
if (SPHINX_FOUND)
   configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in
        ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)
   file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/html)
   add_custom_target(python_doc
        COMMAND ${SPHINX_EXECUTABLE}
                -b html
                -d build/doctrees
                -c .
                -E   # don't use cached files
                ${CMAKE_CURRENT_SOURCE_DIR}
                build/html
        COMMAND ${PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/postprocess_html.py
                ${CMAKE_CURRENT_BINARY_DIR}/build/html
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generate Python API documentation with Sphinx" VERBATIM)
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/html
            DESTINATION bempp/python/doc)
else (SPHINX_FOUND)
    message(WARNING "Sphinx not found; Python documentation will not be created")
endif (SPHINX_FOUND)
