find_package(Git)

# Set version variable
# If git is found and this is a git workdir, then figure out build id
# Stores results in ${PROJECT_NAME}_VERSION and ${PROJECT_NAME}_GITREF
function(set_version VERSION)
  set(${PROJECT_NAME}_VERSION ${VERSION} PARENT_SCOPE)
endfunction()

# Set version variable
# If git is found and this is a git workdir, then figure out build id
# Stores results in ${PROJECT_NAME}_VERSION and ${PROJECT_NAME}_GITREF
function(get_gitref)

  # Tries and get git hash first
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      RESULT_VARIABLE HASH_RESULT
      OUTPUT_VARIABLE GIT_HASH
      ERROR_QUIET
    )
    if(HASH_RESULT EQUAL 0)
      string(STRIP "${GIT_HASH}" GIT_HASH)
    else()
      set(GIT_HASH "NA")
    endif()
  else()
    set(GIT_HASH "NA")
  endif()

  set(${PROJECT_NAME}_GITREF ${GIT_HASH} PARENT_SCOPE)

endfunction()
