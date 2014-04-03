# This file patches two headers from dune-localfunction.
# We use regex rather than patches so that the patch work for subsequent versions.

function(apply_regex FILENAME OLD_PATTERN NEW_PATTERN)
  file(READ ${FILENAME} data)
  string(REGEX REPLACE ${OLD_PATTERN} ${NEW_PATTERN} data "${data}")
  file(WRITE ${FILENAME} "${data}")
endfunction()

if(NOT ROOT)
  message(FATAL_ERROR "ROOT variable not defined")
endif()

# The patches cast 1.0 to an actual type.
apply_regex(
  ${ROOT}/dune/localfunctions/lagrange/pk2d/pk2dlocalbasis.hh
  "product( *)\\*=( *-?)1\\.0"
  "product\\1*=\\2D(1)"
)

# The patches cast 1.0 to an actual type.
apply_regex(
  ${ROOT}/dune/localfunctions/raviartthomas/raviartthomas02d/raviartthomas02dlocalbasis.hh
  "in\\[(0|1)\\]( *)-( *)1\\.0"
  "in[\\1]\\2-\\3D(1)"
)
