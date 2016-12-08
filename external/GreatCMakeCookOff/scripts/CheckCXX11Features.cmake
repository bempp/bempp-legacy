# Checks for C++11 features
#
# USAGE: There are two functions
#
# cxx11_find_all_features(OUTPUT_VARIABLE)
# This function returns a variable with all possible features.
#
# cxx11_feature_check([feature feature] [REQUIRED [feature feature]])
# If no arguments are provided, then checks all available features
# Features appeacing before REQUIRED are optional.
# If arguments are provided and those features are available, sets
# the variable HAS_CXX11_FEATURENAME, where FEATURENAME is the input in capital letters.
# Fails if required feature are not available
#
# For possible features, please print out the result from the first function.
#
# Original script by Rolf Eike Beer
# Modifications by Andreas Weis
# Further Modifications by RSDT@UCL
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)

set(CPP11_FEATURE_CHECK_DIR ${CMAKE_CURRENT_LIST_DIR}/cpp11
    CACHE INTERNAL "c++11 file directory")

MACRO(cxx11_check_single_feature FEATURE_NAME FEATURE_NUMBER RESULT_VAR)
	IF (NOT DEFINED ${RESULT_VAR})
    SET(_bindir "${CMAKE_BINARY_DIR}/cxx11_feature_tests/cxx11_${FEATURE_NAME}")

		IF (${FEATURE_NUMBER})
      SET(_SRCFILE_BASE ${CPP11_FEATURE_CHECK_DIR}/${FEATURE_NAME}-N${FEATURE_NUMBER})
			SET(_LOG_NAME "\"${FEATURE_NAME}\" (N${FEATURE_NUMBER})")
		ELSE (${FEATURE_NUMBER})
      SET(_SRCFILE_BASE ${CPP11_FEATURE_CHECK_DIR}/${FEATURE_NAME})
			SET(_LOG_NAME "\"${FEATURE_NAME}\"")
		ENDIF (${FEATURE_NUMBER})
		MESSAGE(STATUS "Checking C++11 support for ${_LOG_NAME}")

		SET(_SRCFILE "${_SRCFILE_BASE}.cpp")
		SET(_SRCFILE_FAIL "${_SRCFILE_BASE}_fail.cpp")
		SET(_SRCFILE_FAIL_COMPILE "${_SRCFILE_BASE}_fail_compile.cpp")

		IF (CROSS_COMPILING)
			try_compile(${RESULT_VAR} "${_bindir}" "${_SRCFILE}")
			IF (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL})
				try_compile(${RESULT_VAR} "${_bindir}_fail" "${_SRCFILE_FAIL}")
			ENDIF (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL})
		ELSE (CROSS_COMPILING)
			try_run(_RUN_RESULT_VAR _COMPILE_RESULT_VAR
					"${_bindir}" "${_SRCFILE}")
			IF (_COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR)
				SET(${RESULT_VAR} TRUE)
			ELSE (_COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR)
				SET(${RESULT_VAR} FALSE)
			ENDIF (_COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR)
			IF (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL})
				try_run(_RUN_RESULT_VAR _COMPILE_RESULT_VAR
						"${_bindir}_fail" "${_SRCFILE_FAIL}")
				IF (_COMPILE_RESULT_VAR AND _RUN_RESULT_VAR)
					SET(${RESULT_VAR} TRUE)
				ELSE (_COMPILE_RESULT_VAR AND _RUN_RESULT_VAR)
					SET(${RESULT_VAR} FALSE)
				ENDIF (_COMPILE_RESULT_VAR AND _RUN_RESULT_VAR)
			ENDIF (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL})
		ENDIF (CROSS_COMPILING)
		IF (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL_COMPILE})
			try_compile(_TMP_RESULT "${_bindir}_fail_compile" "${_SRCFILE_FAIL_COMPILE}")
			IF (_TMP_RESULT)
				SET(${RESULT_VAR} FALSE)
			ELSE (_TMP_RESULT)
				SET(${RESULT_VAR} TRUE)
			ENDIF (_TMP_RESULT)
		ENDIF (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL_COMPILE})

		IF (${RESULT_VAR})
			MESSAGE(STATUS "Checking C++11 support for ${_LOG_NAME} -- works")
		ELSE (${RESULT_VAR})
			MESSAGE(STATUS "Checking C++11 support for ${_LOG_NAME} -- not supported")
		ENDIF (${RESULT_VAR})
		SET(${RESULT_VAR} ${${RESULT_VAR}} CACHE INTERNAL "C++11 support for ${_LOG_NAME}")
	ENDIF (NOT DEFINED ${RESULT_VAR})
ENDMACRO(cxx11_check_single_feature)

# Find list of all features
function(cxx11_find_all_features outvar)
  FILE(GLOB ALL_CPP11_FEATURE_FILES "${CPP11_FEATURE_CHECK_DIR}/*.cpp")
  set(OUTPUT_VARIABLES)
  foreach(filename ${ALL_CPP11_FEATURE_FILES})
    get_filename_component(filename ${filename} NAME_WE)
    string(REGEX REPLACE "_fail_compile" "" filename "${filename}")
    string(REGEX REPLACE "_fail" "" filename "${filename}")
    string(REGEX REPLACE "-N[0-9]*" "" filename "${filename}")
    set(OUTPUT_VARIABLES ${OUTPUT_VARIABLES} ${filename})
  endforeach()
  FILE(GLOB obsolete "${CPP11_FEATURE_CHECK_DIR}/obsolete/*.cpp")
  foreach(filename ${obsolete})
    get_filename_component(filename ${filename} NAME_WE)
    string(REGEX REPLACE "_fail_compile" "" filename "${filename}")
    string(REGEX REPLACE "_fail" "" filename "${filename}")
    string(REGEX REPLACE "-N[0-9]*" "" filename "${filename}")
    set(OUTPUT_VARIABLES ${OUTPUT_VARIABLES} "obsolete/${filename}")
  endforeach()
  list(REMOVE_DUPLICATES OUTPUT_VARIABLES)
  set(${outvar} ${OUTPUT_VARIABLES} PARENT_SCOPE)
endfunction()

# Parses input and separates into arguments before REQUIRED and after REQUIRED.
# Arguments before REQUIRED are OPTIONALS.
# Arguments after REQUIRED are REQUIRED.
# If no arguments, then sets output OPTIONALS to ALLFEATURES.
function(parse_input_features ALLFEATURES OPTIONALS REQUIRED ERRORS)

  if("${ARGN}" STREQUAL "")
    set(${OPTIONALS} ${ALLFEATURES} PARENT_SCOPE)
    set(${REQUIRED} "" PARENT_SCOPE)
  else()
    set(REQUIRED_FEATURES)
    set(OPTIONAL_FEATURES)
    set(UNKNOWN_FEATURES)
    set(result_type OPTIONAL_FEATURES)
    foreach(feature ${ARGN})
      string(COMPARE EQUAL "${feature}" "REQUIRED" avoid_cmake_warning)
      if(avoid_cmake_warning)
        set(result_type REQUIRED_FEATURES)
      else()
        list(FIND ALLFEATURES ${feature} feature_was_found)

        if(feature_was_found EQUAL -1)
          list(APPEND UNKNOWN_FEATURES ${feature})
        else()
          list(APPEND ${result_type} ${feature})
        endif()
      endif()
    endforeach()

    set(${OPTIONALS} ${OPTIONAL_FEATURES} PARENT_SCOPE)
    set(${REQUIRED} ${REQUIRED_FEATURES} PARENT_SCOPE)
    set(${ERRORS} ${UNKNOWN_FEATURES} PARENT_SCOPE)
  endif("${ARGN}" STREQUAL "")
endfunction(parse_input_features)

# Figures out name and number of feature
# then calls macro that does the work
macro(_figure_out_cxx11_feature current_feature)
  # Find set of files that match current_feature, excepting _fail and _fail_compile.
  file(GLOB ALL_FEATURE_FILES "${CPP11_FEATURE_CHECK_DIR}/${current_feature}*.cpp")
  foreach(filename ${ALL_FEATURE_FILES})
    if(filename MATCHES "_fail")
      list(REMOVE_ITEM ALL_FEATURE_FILES ${filename})
    endif()
  endforeach()

  list(LENGTH ALL_FEATURE_FILES NFILES)
  if(NOT ${NFILES} EQUAL 1)
    message(FATAL_ERROR "[c++11] Expected to find only one feature. Found ${NFILES} -- ${ALL_FEATURE_FILES}.")
  endif(NOT ${NFILES} EQUAL 1)

  # Now we know which file corresponds to option.
  get_filename_component(basename ${ALL_FEATURE_FILES} NAME_WE)
  # If has feature number, extract it
  set(number "")
  if(basename MATCHES "-N[0-9]*$")
    string(REGEX REPLACE "${current_feature}-N" "" number "${basename}")
  endif()
  # Then call macro
  string(TOUPPER ${current_feature} UPPER_OPTIONAL)
  string(REGEX REPLACE "/" "_" VARNAME "HAS_CXX11_${UPPER_OPTIONAL}")
  cxx11_check_single_feature(${current_feature} "${number}" ${VARNAME})
endmacro(_figure_out_cxx11_feature)

function(cxx11_feature_check)

  # find all features
  cxx11_find_all_features(ALL_CPP11_FEATURES)

  # Parses input to this function.
  parse_input_features("${ALL_CPP11_FEATURES}" OPTIONALS REQUIRED ERRORS ${ARGN})
  if(NOT ${ERRORS} STREQUAL "")
    message(STATUS "[c++11] The following features are unknown: ${ERRORS}.")
  endif()

  # MinGW has not implemented std::random_device fully yet. Unfortunately, this can only be detected
  # by running a program which tries to call std::random_device. However that generates an error that
  # is *not* caught by CMake's try_run.
  if(MSYS)
    list(REMOVE_ITEM OPTIONALS "random_device")
    list(FIND REQUIRED "random_device" feature_was_found)
    if(NOT feature_was_found EQUAL "-1")
      message(FATAL_ERROR "[c++1] MSYS does not implement Random devices fully.\n"
                          "       It cannot be required on this system.")
    endif()
  endif()

  # Check optional features
  foreach(current_feature ${OPTIONALS})
    _figure_out_cxx11_feature(${current_feature})
  endforeach(current_feature ${ARGN})

  # Check required features
  foreach(current_feature ${REQUIRED})
    _figure_out_cxx11_feature(${current_feature})
    set(VARNAME HAS_CXX11_${UPPER_OPTIONAL})
    if(NOT ${VARNAME})
      message(FATAL_ERROR "[c++11] Required feature ${current_feature} is not available.")
    endif(NOT ${VARNAME})
  endforeach(current_feature ${REQUIRED})

endfunction(cxx11_feature_check)
