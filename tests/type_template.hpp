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

#ifndef bempp_type_template_hpp
#define bempp_type_template_hpp

#include "bempp/common/config_data_types.hpp"

#include <boost/mpl/list.hpp>
#include <complex>

#ifdef ENABLE_SINGLE_PRECISION
#  ifdef ENABLE_DOUBLE_PRECISION
typedef boost::mpl::list<float, double, std::complex<float>, std::complex<double> >
numeric_types;
typedef boost::mpl::list<float, double>
real_numeric_types;
typedef boost::mpl::list<std::complex<float>, std::complex<double> >
complex_numeric_types;
#  else
typedef boost::mpl::list<float, std::complex<float> >
numeric_types;
typedef boost::mpl::list<float>
real_numeric_types;
typedef boost::mpl::list<std::complex<float> >
complex_numeric_types;
#  endif
#else
#  ifdef ENABLE_DOUBLE_PRECISION
typedef boost::mpl::list<double, std::complex<double> >
numeric_types;
typedef boost::mpl::list<double>
real_numeric_types;
typedef boost::mpl::list<std::complex<double> >
complex_numeric_types;
#  else
typedef boost::mpl::list<>
numeric_types;
typedef boost::mpl::list<>
real_numeric_types;
typedef boost::mpl::list<>
complex_numeric_types;
#  endif
#endif

#ifdef ENABLE_COMPLEX_BASIS_FUNCTIONS
typedef numeric_types basis_function_types;
#else
typedef real_numeric_types basis_function_types;
#endif

#ifdef ENABLE_COMPLEX_KERNELS
typedef numeric_types kernel_types;
#else
typedef real_numeric_types kernel_types;
#endif

#if defined(ENABLE_COMPLEX_KERNELS) || defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
typedef numeric_types result_types;
typedef complex_numeric_types complex_result_types;
#else
typedef real_numeric_types result_types;
typedef boost::mpl::list<> complex_result_types;
#endif

#endif
