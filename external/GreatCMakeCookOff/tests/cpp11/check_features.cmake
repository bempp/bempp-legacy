# This is a C++ feature
enable_language(CXX)

# Include check feature.
include(${cookoff_path}/scripts/CheckCXX11Features.cmake)

# Now call function to test.
cxx11_feature_check("auto")
cxx11_feature_check()
