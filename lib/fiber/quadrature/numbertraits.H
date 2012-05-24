// Copyright (C) 2009-2010 Matthias Messner, Michael Messner, Franz
// Rammerstorfer, Peter Urthaler
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

// Note: This file is originally part of the HyENA project
// (http://portal.tugraz.at/portal/page/portal/Files/i2610/files/Forschung/Software/HyENA/html/index.html)
// and has been relicensed with permission by the HyENA authors. This does not affect the license of any
// other part of HyENA.

/**
 * @file   numbertraits.H
 * @ingroup traits
 *
 * The @p NumberTraits - A structure that provides traits and member
 * functions for type @p T that make it possible to write templates that 
 * work on real and complex numbers.  Mostly used by linear algebra classes.
 * Feel free to add your type!  
 *
 * @tparam T type
 *
 * @author Rf
 * @date   created:     02.10.09
 *         last change: 02.11.09
 */
#ifndef numbertraits_h
#define numbertraits_h



// system includes
#include <complex>
#include <cmath>



// base struct - not defined
template<typename T> struct NumberTraits;



// specialisation for int
template<>
struct NumberTraits<int>
{
	// some typedefs
	typedef int            value_type;
	typedef value_type      base_type;
	typedef value_type  argument_type;
	typedef double         float_type;

	// flag for complex numbers.
	static const bool is_complex = false;

	// return zero.
  static inline value_type zero() { return 0; }

	// return one.
  static inline value_type one() { return 1; }

	// return tolerance (always real value).
  static inline base_type tolerance() { return 0; }

	// return the real part of the given number.
	static inline	base_type real(argument_type x) { return x; }
	
	// return the imaginary part of the given number.
	static inline	base_type imag(argument_type x) { return 0; }
	
	// Return the complex-conjugate of the given number.
	static inline	value_type conj(argument_type x) { return x; }
	
	// return the absolute value of the given number.
	static inline	base_type abs(argument_type x) { return std::fabs(x); }
};


// specialisation for int
template<>
struct NumberTraits<unsigned>
{
	// some typedefs
	typedef unsigned       value_type;
	typedef value_type      base_type;
	typedef value_type  argument_type;
	typedef double         float_type;

	// flag for complex numbers.
	static const bool is_complex = false;

	// return zero.
  static inline value_type zero() { return 0; }

	// return one.
  static inline value_type one() { return 1; }

	// return tolerance (always real value).
  static inline base_type tolerance() { return 0; }

	// return the real part of the given number.
	static inline	base_type real(argument_type x) { return x; }
	
	// return the imaginary part of the given number.
	static inline	base_type imag(argument_type x) { return 0; }
	
	// Return the complex-conjugate of the given number.
	static inline	value_type conj(argument_type x) { return x; }
	
	// return the absolute value of the given number.
	static inline	base_type abs(argument_type x) { return std::fabs(x); }
};

// specialisation for float	
template<>
struct NumberTraits<float>
{
	// some typedefs
	typedef float          value_type;
	typedef value_type      base_type;
	typedef value_type  argument_type;
	typedef value_type     float_type;

	// flag for complex numbers.
	static const bool is_complex = false;

	// return zero.
  static inline value_type zero() { return 0.0; }

	// return one.
  static inline value_type one() { return 1.0; }

	// return tolerance (always real value).
  static inline base_type tolerance() { return 1.e-3; }

	// return the real part of the given number.
	static inline	base_type real(argument_type x) { return x; }
	
	// return the imaginary part of the given number.
	static inline	base_type imag(argument_type x) { return 0.0; }
	
	// Return the complex-conjugate of the given number.
	static inline	value_type conj(argument_type x) { return x; }
	
	// return the absolute value of the given number.
	static inline	base_type abs(argument_type x) { return std::fabs(x); }
};



// specialisation for double
template<>
struct NumberTraits<double>
{
	// some typedefs
	typedef double          value_type;
	typedef value_type       base_type;
	typedef value_type   argument_type;
	typedef value_type      float_type;

	// flag for complex numbers.
	static const bool is_complex = false;

	// return zero.
  static inline value_type zero() { return 0.0; }

	// return one.
  static inline value_type one() { return 1.0; }

	// return tolerance (always real value).
  static inline base_type tolerance() { return 1.e-6; }

	// return the real part of the given number.
	static inline	base_type real(argument_type x) { return x; }
	
	// return the imaginary part of the given number.
	static inline	base_type imag(argument_type x) { return 0.0; }
	
	// Return the complex-conjugate of the given number.
	static inline	value_type conj(argument_type x) { return x; }
	
	// return the absolute value of the given number.
	static inline	base_type abs(argument_type x) { return std::fabs(x); }
};



// specialisation for complex<float>
template<>
struct NumberTraits<std::complex<float> >
{
	// some typedefs
	typedef std::complex<float>  value_type;
	typedef float                 base_type;
	typedef const value_type& argument_type;
	typedef value_type           float_type;

	// flag for complex numbers.
	static const bool is_complex = true;

	// return zero.
  static inline value_type zero() { return value_type(0.0,0.0); }

	// return one.
  static inline value_type one() { return value_type(1.0,0.0); }

	// return tolerance (always real value).
  static inline base_type tolerance() { return 1.e-3; }

	// return the real part of the given number.
	static inline	base_type real(argument_type x) { return std::real(x); }
	
	// return the imaginary part of the given number.
	static inline	base_type imag(argument_type x) { return std::imag(x); }
	
	// Return the complex-conjugate of the given number.
	static inline	value_type conj(argument_type x) { return std::conj(x); }
	
	// return the absolute value of the given number.
	static inline	base_type abs(argument_type x) { return std::abs(x); }
};



// specialisation for complex<double>
template <>
struct NumberTraits<std::complex<double> >
{
	// some typedefs
	typedef std::complex<double>  value_type;
	typedef double                 base_type;
	typedef const value_type&  argument_type;
	typedef value_type            float_type;

	// flag for complex numbers.
	static const bool is_complex = true;

	// return zero.
  static inline value_type zero() { return value_type(0.0,0.0); }

	// return one.
  static inline value_type one() { return value_type(1.0,0.0); }

	// return tolerance (always real value).
  static inline base_type tolerance() { return  1.e-6; }

	// return the real part of the given number.
	static inline	base_type real(argument_type x) { return std::real(x); }
	
	// return the imaginary part of the given number.
	static inline	base_type imag(argument_type x) { return std::imag(x); }
	
	// Return the complex-conjugate of the given number.
	static inline	value_type conj(argument_type x) { return std::conj(x); }
	
	// return the absolute value of the given number.
	static inline	base_type abs(argument_type x) { return std::abs(x); }
};



#endif

