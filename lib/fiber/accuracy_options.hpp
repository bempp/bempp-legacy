#ifndef fiber_accuracy_options_hpp
#define fiber_accuracy_options_hpp

#include "quadrature_options.hpp"

namespace Fiber
{

/** \brief Options controlling quadrature accuracy. */
struct AccuracyOptions
{
public:
    QuadratureOptions regular;
    QuadratureOptions singular;
};

} // namespace Fiber

#endif
