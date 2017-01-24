if(NOT Julia_FOUND)
  find_program(Julia_EXECUTABLE julia DOC "Julia executable")
  if(Julia_EXECUTABLE)
      execute_process(
          COMMAND ${Julia_EXECUTABLE} -E "Pkg.dir()"
          OUTPUT_VARIABLE Julia_INSTALL_DIR
          RESULT_VARIABLE RESULT
      )
      if(RESULT EQUAL 0)
          string(REGEX REPLACE "\"" "" Julia_INSTALL_DIR ${Julia_INSTALL_DIR})
          set(Julia_INSTALL_DIR ${Julia_INSTALL_DIR}
              CACHE PATH "Path where Julia packages should be installed.")
      endif()
      execute_process(
          COMMAND ${Julia_EXECUTABLE} --version
          OUTPUT_VARIABLE Julia_VERSION_STRING
          RESULT_VARIABLE RESULT
      )
      if(RESULT EQUAL 0)
        string(REGEX REPLACE ".*([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1"
            Julia_VERSION_STRING ${Julia_VERSION_STRING})
      endif()
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Julia
      REQUIRED_VARS Julia_EXECUTABLE Julia_INSTALL_DIR
      VERSION_VAR Julia_VERSION_STRING
  )
endif()
