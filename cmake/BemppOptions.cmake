# Options (can be modified by user)
option(WITH_TESTS "Compile unit tests (can be run with 'make test')" ON)
option(WITH_AHMED "Link to the AHMED library to enable ACA mode assembly)" OFF)
option(WITH_TRILINOS "Use sparse linear solvers from the Trilinos library)" OFF)

set(DUNE_PATH ${CMAKE_SOURCE_DIR}/contrib/dune CACHE PATH "Path to Dune")
