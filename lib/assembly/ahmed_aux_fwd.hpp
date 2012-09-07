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


#ifndef bempp_ahmed_aux_fwd_hpp
#define bempp_ahmed_aux_fwd_hpp

#include "../common/common.hpp"

#include "bempp/common/config_ahmed.hpp"

/** \cond FORWARD_DECL */
template <class T1, class T2> class bemblcluster;
template <class T> class mblock;
class blcluster;
/** \endcond */

namespace std
{

/** \cond FORWARD_DECL */
template <typename _Tp> class complex;
/** \endcond */

} // namespace std

#ifdef AHMED_USES_STD_COMPLEX
typedef std::complex<float> scomp;
typedef std::complex<double> dcomp;
#else
template <typename T> class comp;
typedef comp<float> scomp;
typedef comp<double> dcomp;
#endif

namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename CoordinateType> struct AhmedDofWrapper;
template <typename T> class ExtendedBemCluster;
/** \endcond */

// Casts.

// If Ahmed is compiled with STD_CMPLX undefined, it uses a
// nonstandard complex type. To all probability, its binary
// representation is the same as that of the complex type provided by
// STL.

template <typename T>
struct AhmedTypeTraits
{
    typedef T Type;
};

template <>
struct AhmedTypeTraits<std::complex<float> >
{
    typedef scomp Type;
};

template <>
struct AhmedTypeTraits<std::complex<double> >
{
    typedef dcomp Type;
};

template <typename T>
inline typename AhmedTypeTraits<T>::Type* ahmedCast(T* x) {
#ifdef AHMED_USES_STD_COMPLEX
    return x;
#else
    return reinterpret_cast<typename AhmedTypeTraits<T>::Type*>(x);
#endif 
}

float ahmedCast(float x);
double ahmedCast(double x);
scomp ahmedCast(std::complex<float> x);
dcomp ahmedCast(std::complex<double> x);

} // namespace Bempp

#endif
