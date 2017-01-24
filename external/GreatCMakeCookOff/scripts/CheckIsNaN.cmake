# Tries and find which isnan to use
# defines ISNAN_VARIATION, which should be included in a configuration files somewher as
#   @ISNAN_VARIATION@
if(ISNAN_VARIATION)
  return()
endif(ISNAN_VARIATION)

include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES(
  "#include <cmath>\nint main() { bool a = std::isnan(0e0); return 0; }\n" 
  CXX_HAS_STD_ISNAN)

if(NOT CXX_HAS_STD_ISNAN)
  CHECK_CXX_SOURCE_COMPILES(
    "#include <math.h>\nint main() { bool a = isnan(0e0); return 0; }\n" 
    CXX_HAS_ISNAN)
endif(NOT CXX_HAS_STD_ISNAN)

if(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN)
  CHECK_CXX_SOURCE_COMPILES(
    "#include <math.h>\nint main() { bool a = _isnan(0e0); return 0; }\n" 
    CXX_HAS___ISNAN)
endif(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN)

if(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN)
  CHECK_CXX_SOURCE_COMPILES(
    "# include <float.h>\nint main() { bool a = _isnan(0e0); return 0; }\n" 
    CXX_HAS_FLOAT_H_ISNAN)
endif(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN)

if(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN AND NOT CXX_HAS___ISNAN AND NOT CXX_HAS_FLOAT_H_ISNAN)
  message(FATAL_ERROR "[isnan] could not find standard function on this OS.")
endif(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN AND NOT CXX_HAS___ISNAN AND NOT CXX_HAS_FLOAT_H_ISNAN)

if(CXX_HAS_STD_ISNAN)
  set(ISNAN_HEADERS "#include <cmath>")
  set(ISNAN_VARIATION "std::isnan")
elseif(CXX_HAS_ISNAN)
  set(ISNAN_HEADERS "#include <math.h>")
  set(ISNAN_VARIATION "::isnan")
elseif(CXX_HAS___ISNAN)
  set(ISNAN_HEADERS "#include <math.h>")
  set(ISNAN_VARIATION "__isnan")
elseif(CXX_HAS_FLOAT_H_ISNAN)
  set(ISNAN_HEADERS "#include <float.h>")
  set(ISNAN_VARIATION "_isnan")
else()
  message(FATAL_ERROR "AM HERE")
endif()
if(ISNAN_VARIATION)
  set(ISNAN_VARIATION ${ISNAN_VARIATION} CACHE INTERNAL "Definition for isnan\n")
  set(ISNAN_HEADERS ${ISNAN_HEADERS} CACHE INTERNAL "Headers containing isnan definition\n")
endif(ISNAN_VARIATION)
