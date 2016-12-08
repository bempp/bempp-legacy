find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()

set(EXTERNAL_ROOT "${PROJECT_BINARY_DIR}/../lookup_eigen/external")
include(PackageLookup)
lookup_package(Eigen3 DOWNLOAD_BY_DEFAULT CHECK_EXTERNAL)

if(NOT Eigen3_BUILT_AS_EXTERNAL_PROJECT)
  message(FATAL_ERROR "not an external project")
endif()
