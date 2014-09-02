# Examples
add_executable(tutorial_dirichlet tutorial_dirichlet.cpp)
target_link_libraries(tutorial_dirichlet libbempp)

install(TARGETS tutorial_dirichlet
    EXPORT BemppTargets
    RUNTIME
    DESTINATION ${RUNTIME_INSTALL_PATH}/bempp/examples)

install(FILES tutorial_dirichlet.cpp DESTINATION ${SHARE_INSTALL_PATH}/bempp/examples/cpp)

