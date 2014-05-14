# Checks for the ALUGrid library: http://aam.mathematik.uni-freiburg.de/IAM/Research/alugrid/
# Sets the following variables
# - ALUGrid_LIBRARIES
# - ALUGrid_INCLUDE_DIR
find_library(ALUGrid_LIBRARIES NAMES alugrid)
find_path(ALUGrid_INCLUDE_DIR alugrid_2d.h
    PATH_SUFFIXES ALUGrid alugrid
)
find_path(ALUGrid_GITTER_INCLUDE_DIR gitter_impl.h
    PATH_SUFFIXES ALUGrid/serial alugrid/serial
)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
  ALUGrid
  DEFAULT_MSG
  ALUGrid_LIBRARIES ALUGrid_INCLUDE_DIR
  ALUGrid_GITTER_INCLUDE_DIR
)
