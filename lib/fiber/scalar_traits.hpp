// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef fiber_scalar_traits_hpp
#define fiber_scalar_traits_hpp

#include "../common/common.hpp"

#include <complex>

namespace Fiber
{

template <typename T>
struct ScalarTraits
{
    // If you get a compilation error here, you are probably trying to use an
    // unsupported floating-point type. The supported types are: float, double,
    // std::complex<float> and std::complex<double>.
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
    // If you get a compilation error here, chances are that you are trying to
    // mix floating-point numbers with different precisions (e.g. float and
    // double or std::complex<float> and double). BEM++ does not support this.
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
