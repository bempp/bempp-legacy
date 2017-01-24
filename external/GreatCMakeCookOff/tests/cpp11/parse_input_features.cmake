# This is a C++ feature
enable_language(CXX)

# Include check feature.
include(${cookoff_path}/scripts/CheckCXX11Features.cmake)

# Now call function with fake input
set(ALL_FEATURES "first;second;third")

# No input OPTIONALS == ALL_FEATURES
parse_input_features("${ALL_FEATURES}" OPTIONALS REQUIRED ERRORS "")
if(NOT "${OPTIONALS}" STREQUAL "${ALL_FEATURES}")
  message(FATAL_ERROR "OPTIONALS results should contain everything.") 
endif(NOT "${OPTIONALS}" STREQUAL "${ALL_FEATURES}")
if(NOT "${REQUIRED}" STREQUAL "")
  message(FATAL_ERROR "REQUIRED should be empty.") 
endif(NOT "${REQUIRED}" STREQUAL "")
if(NOT "${ERRORS}" STREQUAL "")
  message(FATAL_ERROR "ERRORS should be empty.") 
endif(NOT "${ERRORS}" STREQUAL "")

# Single input without REQUIRED
parse_input_features("${ALL_FEATURES}" OPTIONALS REQUIRED ERRORS "second")
if(NOT "${OPTIONALS}" STREQUAL "second")
  message(FATAL_ERROR "OPTIONALS results should second.") 
endif()
if(NOT "${REQUIRED}" STREQUAL "")
  message(FATAL_ERROR "REQUIRED should be empty.") 
endif()
if(NOT "${ERRORS}" STREQUAL "")
  message(FATAL_ERROR "ERRORS should be empty.") 
endif()

# Single error input without REQUIRED
parse_input_features("${ALL_FEATURES}" OPTIONALS REQUIRED ERRORS "none")
if(NOT "${OPTIONALS}" STREQUAL "")
  message(FATAL_ERROR "OPTIONALS results should be empty.") 
endif()
if(NOT "${REQUIRED}" STREQUAL "")
  message(FATAL_ERROR "REQUIRED should be empty.") 
endif()
if(NOT "${ERRORS}" STREQUAL "none")
  message(FATAL_ERROR "ERRORS should be none.") 
endif()

# Single valid input with REQUIRED
parse_input_features("${ALL_FEATURES}" OPTIONALS REQUIRED ERRORS REQUIRED first)
if(NOT "${OPTIONALS}" STREQUAL "")
  message(FATAL_ERROR "OPTIONALS results should be empty.") 
endif()
if(NOT "${REQUIRED}" STREQUAL "first")
  message(FATAL_ERROR "REQUIRED should be first.") 
endif()
if(NOT "${ERRORS}" STREQUAL "")
  message(FATAL_ERROR "ERRORS should be empty.") 
endif()

# one of each
parse_input_features("${ALL_FEATURES}" OPTIONALS REQUIRED ERRORS second third REQUIRED first none)
if(NOT "${OPTIONALS}" STREQUAL "second;third")
  message(FATAL_ERROR "OPTIONALS results should be second;optional.") 
endif()
if(NOT "${REQUIRED}" STREQUAL "first")
  message(FATAL_ERROR "REQUIRED should be first.") 
endif()
if(NOT "${ERRORS}" STREQUAL "none")
  message(FATAL_ERROR "ERRORS should be none.") 
endif()
