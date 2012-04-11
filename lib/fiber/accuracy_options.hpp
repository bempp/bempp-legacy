#ifndef fiber_accuracy_options_hpp
#define fiber_accuracy_options_hpp

#include "quadrature_options.hpp"

namespace Fiber
{

/** \brief Options controlling quadrature accuracy. */
struct AccuracyOptions
{
public:
    /** \brief Integration of regular functions on pairs of elements */
    QuadratureOptions doubleRegular;
    /** \brief Integration of singular functions on pairs of elements */
    QuadratureOptions doubleSingular;
};

} // namespace Fiber

#endif
