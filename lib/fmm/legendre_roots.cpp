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

#include "legendre_roots.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include <boost/math/special_functions/legendre.hpp>
#include <iostream>		// std::cout

namespace Bempp
{

// P'_n(x): Derivative of the n_{th} order Legendre polynomial w.r.t. x
template <typename ValueType> // must be a real type
ValueType diff_legendre_p(int n, ValueType x)
{
	using namespace boost::math;
	return n*( x*legendre_p(n, x) - legendre_p(n-1, x) ) / (x*x - 1);
}

// from http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
// find the Gauss Legendre quadrature points, P_n(xi) = 0
template <typename ValueType>  // must be a real type
LegendreRoots<ValueType>::LegendreRoots(unsigned int N, ValueType *roots, ValueType *weights)
{
	using namespace boost::math;
	ValueType pi = boost::math::constants::pi<ValueType>();

	roots[int(N/2)] = 0; // take advantage of symmetry
	for (unsigned int i = 1; i <= int(N/2); i++) {
		ValueType x, x1;
		x = cos(pi * (i - .25) / (N + .5)); // guess root position
		int iter=100;
		do { // apply Newton-Raphson method to find roots
			x1 = x;
			x -= legendre_p(N, x) / diff_legendre_p(N, x);
		} while (x != x1 && --iter); // well-behaved function, convergence guaranteed
		roots[i-1] =  x;
		roots[N-i] = -x;
 	}

	for (unsigned int i = 1; i <= int(N/2)+1; i++) {
		ValueType x = roots[i-1];
		ValueType diffPx = diff_legendre_p(N, x);
		weights[i-1] = 2 / ((1 - x*x) * diffPx*diffPx);
		weights[N-i] = weights[i-1];
	}

/*	for (unsigned int i = 1; i <= N; i++)
		std::cout << roots[i - 1] << ", ";, unsigned int levels
	std::cout <<std::endl;
 
	for (unsigned int i = 1; i <= N; i++)
		std::cout << weights[i - 1] << ", ";
	std::cout <<std::endl;
*/
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_REAL_ONLY(
    LegendreRoots);

} // namespace Bempp
