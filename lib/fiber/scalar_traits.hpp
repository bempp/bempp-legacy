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
    typedef std::complex<float> ComplexType;
};

template <>
struct ScalarTraits<double>
{
    typedef double RealType;
    typedef std::complex<double> ComplexType;
};

template <>
struct ScalarTraits<std::complex<float> >
{
    typedef float RealType;
    typedef std::complex<float> ComplexType;
};

template <>
struct ScalarTraits<std::complex<double> >
{
    typedef double RealType;
    typedef std::complex<double> ComplexType;
};

/** \brief "Larger" of the types U and V. */
template <typename U, typename V>
struct Coercion
{
};

template <>
struct Coercion<float, float>
{
    typedef float Type;
};

template <>
struct Coercion<double, double>
{
    typedef double Type;
};

template <>
struct Coercion<std::complex<float>, std::complex<float> >
{
    typedef std::complex<float> Type;
};

template <>
struct Coercion<std::complex<double>, std::complex<double> >
{
    typedef std::complex<double> Type;
};

template <>
struct Coercion<float, std::complex<float> >
{
    typedef std::complex<float> Type;
};

template <>
struct Coercion<std::complex<float>, float>
{
    typedef std::complex<float> Type;
};

template <>
struct Coercion<double, std::complex<double> >
{
    typedef std::complex<double> Type;
};

template <>
struct Coercion<std::complex<double>, double>
{
    typedef std::complex<double> Type;
};

} // namespace Fiber

#endif
