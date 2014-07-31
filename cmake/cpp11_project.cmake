# Flags to compile the project with C++11
include(AddCPP11Flags)

# Makes it possible to check for C++11 features
include(CheckCXX11Features)

cxx11_feature_check(
    steady_clock obsolete/monotonic_clock
    REQUIRED unique_ptr nullptr
)
if(NOT HAS_CXX11_STEADY_CLOCK AND HAS_CXX11_OBSOLETE_MONOTONIC_CLOCK)
  set(BEMPP_ADD_STEADY_CLOCK_FROM_MONOTONIC_CLOCK TRUE)
endif()
