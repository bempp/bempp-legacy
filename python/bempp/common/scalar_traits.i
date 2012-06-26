%{
#include "common/scalar_traits.hpp"
%}

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

%template() ScalarTraits<float>;
%template() ScalarTraits<double>;
%template() ScalarTraits<std::complex<float> >;
%template() ScalarTraits<std::complex<double> >;

} // namespace Fiber


// we repeat the above for the Bempp namespace, since Swig doesn't always process
// using declarations correctly

namespace Bempp
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

%template() ScalarTraits<float>;
%template() ScalarTraits<double>;
%template() ScalarTraits<std::complex<float> >;
%template() ScalarTraits<std::complex<double> >;

} // namespace Bempp


%{

namespace Fiber
{

template <typename T>
struct PythonScalarTraits
{
};

template <>
struct PythonScalarTraits<float>
{
    enum { numpyType = NPY_FLOAT };
};

template <>
struct PythonScalarTraits<double>
{
    enum { numpyType = NPY_DOUBLE };
};

template <>
struct PythonScalarTraits<std::complex<float> >
{
    enum { numpyType = NPY_CFLOAT };
};

template <>
struct PythonScalarTraits<std::complex<double> >
{
    enum { numpyType = NPY_CDOUBLE };
};

} // namespace Fiber

// this using declaration is OK since it will be only seen by the C++
// compiler, not by SWIG
namespace Bempp
{
using Fiber::PythonScalarTraits;
} // namespace Bempp

%}
