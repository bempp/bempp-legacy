# Try to find Eigen3
# Will define
#
# Eigen3_FOUND
# Eigen3_INCLUDE_DIR

if(NOT Eigen3_INCLUDE_DIR)
    set(Eigen3_DIR ${PROJECT_BINARY_DIR}/external/include)
endif()

find_path(Eigen3_INCLUDE_DIR Eigen/Eigen
    HINTS ${Eigen3_DIR}
    DOC "Eigen3 Directory"
    )

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
  Eigen3 
  DEFAULT_MSG
  Eigen3_INCLUDE_DIR
)
