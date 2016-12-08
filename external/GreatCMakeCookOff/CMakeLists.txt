cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)
project(GreatCMakeCookOff)

option(tests          "Enable testing."                         on)

set(INSTALL_DIRECTORY share/GreatCMakeCookOff )
add_subdirectory(modules)
add_subdirectory(scripts)

install(
  FILES cmake/GreatCMakeCookOffConfig.cmake
  DESTINATION ${INSTALL_DIRECTORY}/cmake COMPONENT dev
)

install(DIRECTORY hunter_recipes DESTINATION ${INSTALL_DIRECTORY})

if(tests)
  enable_testing()
  list(APPEND CMAKE_MODULE_PATH
        ${PROJECT_SOURCE_DIR}/modules ${PROJECT_SOURCE_DIR}/scripts)
  include(TestCMake)
  add_subdirectory(tests)
endif()
