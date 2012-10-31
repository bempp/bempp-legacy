# This script ensures that the AHMED_G77CONVENTION option is properly
# defined even if the user updates manually from 1.0 to 1.0.1 (without
# a fresh "git clone"). If the option is already in cache, this script
# does nothing.

# This script can be removed in 1.1 (if we assume that users will need
# to clone the git repository afresh).

if (NOT DEFINED AHMED_G77CONVENTION)
   file (READ "${CMAKE_BINARY_DIR}/../.options.cfg" OPTIONS_CFG_CONTENTS)
   string(REGEX MATCH "AHMED_with_g77=\\\"[^\n]*\\\"" AHMED_WITH_G77_OPTION ${OPTIONS_CFG_CONTENTS})
   string(REGEX REPLACE ".*\\\"(.*)\\\".*" "\\1" AHMED_WITH_G77_OPTION_VALUE ${AHMED_WITH_G77_OPTION})
   set(AHMED_G77CONVENTION ${AHMED_WITH_G77_OPTION_VALUE} CACHE BOOL "Do BLAS and LAPACK use G77 calling convention?")
endif ()
