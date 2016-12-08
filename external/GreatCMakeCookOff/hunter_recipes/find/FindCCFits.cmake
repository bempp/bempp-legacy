find_path(CCFits_INCLUDE_DIR NAMES CCfits.h PATH_SUFFIXES CCfits)
find_library(CCFits_LIBRARY NAMES CCfits)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
  CCFits  DEFAULT_MSG CCFits_INCLUDE_DIR CCFits_LIBRARY)
