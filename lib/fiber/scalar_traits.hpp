#ifndef fiber_scalar_traits_hpp
#define fiber_scalar_traits_hpp

#include <complex>

namespace Fiber
{

template <typename T>
struct ScalarTraits
{

};

template <>
struct ScalarTraits<float>
{
    typedef float RealType;
};

template <>
struct ScalarTraits<double>
{
    typedef double RealType;
};

template <>
struct ScalarTraits<std::complex<float> >
{
    typedef float RealType;
};

template <>
struct ScalarTraits<std::complex<double> >
{
    typedef double RealType;
};


} // namespace Fiber

#endif
