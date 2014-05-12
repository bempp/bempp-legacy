#ifndef BEMPP_CONFIG_BEMPP_H
#define BEMPP_CONFIG_BEMPP_H
namespace Bempp {
// Makes it possible to use older, pre-release version of C++11 standard
// library, and especially monotonic_clock, rather than steady_clock
// A steady clock is defined from monotonic clock in armadillo_fwd.hpp
#cmakedefine BEMPP_ADD_STEADY_CLOCK_FROM_MONOTONIC_CLOCK
}
#endif
