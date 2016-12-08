file(GLOB scripts *.cmake *.py)
file(GLOB numpy numpy/*.cc numpy/*.c)

install(FILES ${scripts} DESTINATION ${INSTALL_DIRECTORY}/modules)
install(FILES ${numpy} DESTINATION ${INSTALL_DIRECTORY}/modules/numpy)
