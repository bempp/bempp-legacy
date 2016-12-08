# Looks for flags to turn on c99
include(CheckCCompilerFlag)
check_c_compiler_flag(-std=c99 has_std_c99)
if(has_std_c99)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
endif()
