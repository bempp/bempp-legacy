# Add a python_doc target to generate Python API documentation with Sphinx
if (SPHINX_FOUND)
    configure_file(
     ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in
         ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)
    file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/doc/python")
    add_custom_target(python_documentation
        COMMAND ${LOCAL_PYTHON_EXECUTABLE} ${SPHINX_EXECUTABLE}
                -b html
                -d ${CMAKE_CURRENT_BINARY_DIR}/doctrees
                -c ${CMAKE_CURRENT_BINARY_DIR}
                -E   # don't use cached files
                ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_BINARY_DIR}/doc/python/html
        COMMENT "Generate Python API documentation with Sphinx" VERBATIM
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/conf.py
    )
    add_dependencies(documentation python_documentation)
    #    if(EXISTS "${PROJECT_BINARY_DIR}/doc/python")
    #    install(DIRECTORY "${PROJECT_BINARY_DIR}/doc/python"
    #        DESTINATION ${DOC_INSTALL_PATH}/bempp/python)
    #endif()
else (SPHINX_FOUND)
    message(WARNING "Sphinx not found;"
        " Python documentation will not be created")
endif (SPHINX_FOUND)
