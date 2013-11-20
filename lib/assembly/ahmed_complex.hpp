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

#ifndef bempp_ahmed_complex_hpp
#define bempp_ahmed_complex_hpp

#include "../common/common.hpp"
#include <complex>

/** \file
 *
 *  Ahmed's complex numbers */
#include <cmplx.h>

#ifndef AHMED_USES_STD_COMPLEX
inline scomp operator*=(scomp& a, double b)
{
    return operator*=(a, static_cast<float>(b));
}

inline scomp operator/(scomp& a, double b)
{
    return operator/(a, static_cast<float>(b));
}

inline scomp operator/(double a, scomp& b)
{
    return operator/(static_cast<float>(a), b);
}
#endif

namespace Bempp
{

inline float ahmedCast(float x) {
    return x;
}

inline double ahmedCast(double x) {
    return x;
}

inline scomp ahmedCast(std::complex<float> x) {
    return scomp(x.real(), x.imag());
}

inline dcomp ahmedCast(std::complex<double> x) {
    return dcomp(x.real(), x.imag());
}

} // namespace Bempp

#endif
