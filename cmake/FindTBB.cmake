# Finds the thread building block library
# Sets the following variables
# - TBB_LIBRARY
# - TBB_LIBRARY_DEBUG
# - TBB_MALLOC_LIBRARY
# - TBB_MALLOC_LIBRARY_DEBUG
# - TBB_INCLUDE_DIR

find_library(TBB_LIBRARY NAMES libtbb.so libtbb.dylib)
find_library(TBB_LIBRARY_DEBUG NAMES libtbb_debug.so libtbb_debug.dylib)
find_library(TBB_MALLOC_LIBRARY NAMES libtbbmalloc.so libtbbmalloc.dylib)
find_library(TBB_MALLOC_LIBRARY_DEBUG NAMES libtbbmalloc_debug.so libtbbmalloc_debug.dylib)
find_path(TBB_INCLUDE_DIR tbb.h PATH_SUFFIXES tbb)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
  TBB 
  DEFAULT_MSG 
  TBB_LIBRARY TBB_LIBRARY_DEBUG TBB_MALLOC_LIBRARY TBB_MALLOC_LIBRARY_DEBUG
  TBB_INCLUDE_DIR
)
