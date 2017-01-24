find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()

set(EXTERNAL_ROOT "${PROJECT_BINARY_DIR}/../lookup_eigen/external")
include(PackageLookup)
lookup_package(Eigen3)

if(NOT (Eigen3_FOUND OR EIGEN3_FOUND OR Eigen3_BUILT_AS_EXTERNAL_PROJECT))
  message(FATAL_ERROR "not an external project ${Eigen3_BUILT_AS_EXTERNAL_PROJECT} - ${EIGEN3_FOUND} - ${Eigen3_FOUND}")
endif()
