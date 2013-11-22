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

#ifndef bempp_fmm_farfield_function_multiplying_test_hpp
#define bempp_fmm_farfield_function_multiplying_test_hpp

#include "../common/scalar_traits.hpp"
#include "../common/armadillo_fwd.hpp"

namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename ResultType> class FmmTransform;
/** \endcond */

template <typename ResultType>
class FmmFarfieldFunctionMultiplyingTest
{
public:
	typedef ResultType ValueType;
	typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

	FmmFarfieldFunctionMultiplyingTest(const arma::Col<CoordinateType>& khat, 
		const arma::Col<CoordinateType>& nodeCentre,
		const arma::Col<CoordinateType>& nodeSize,
		const FmmTransform<ResultType>& fmmTransform);

	// Number of components of the function's argument
	int argumentDimension() const;
	// Number of components of the function's value
	int resultDimension() const;

	// Evaluate the function at the point "point" and store result in
	// the array "result"
	void evaluate(	const arma::Col<CoordinateType>& point,
				const arma::Col<CoordinateType>& normal,
				arma::Col<ValueType>& result) const;
private: 
	const arma::Col<CoordinateType> &m_khat;
	const arma::Col<CoordinateType> &m_nodeCentre;
	const arma::Col<CoordinateType> &m_nodeSize;
	const FmmTransform<ResultType>& m_fmmTransform;
};

} // namespace Bempp

#endif
