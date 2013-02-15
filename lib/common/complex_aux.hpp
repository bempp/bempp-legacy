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

#ifndef bempp_complex_aux_hpp
#define bempp_complex_aux_hpp

#include "scalar_traits.hpp"

namespace Fiber
{

template <typename T>
inline typename ScalarTraits<T>::RealType realPart(const T& x)
{
    return x;
}

template <>
inline float realPart(const std::complex<float>& x)
{
    return x.real();
}

template <>
inline double realPart(const std::complex<double>& x)
{
    return x.real();
}

template <typename T>
inline typename ScalarTraits<T>::RealType imagPart(const T& x)
{
    return 0.;
}

template <>
inline float imagPart(const std::complex<float>& x)
{
    return x.imag();
}

template <>
inline double imagPart(const std::complex<double>& x)
{
    return x.imag();
}

template <typename T>
inline T conj(const T& x)
{
    return x;
}

template <>
inline std::complex<float> conj(const std::complex<float>& x)
{
    return std::conj(x);
}

template <>
inline std::complex<double> conj(const std::complex<double>& x)
{
    return std::conj(x);
}

} // namespace Fiber

namespace Bempp
{
using Fiber::realPart;
using Fiber::imagPart;
using Fiber::conj;
} // namespace Bempp

#endif
