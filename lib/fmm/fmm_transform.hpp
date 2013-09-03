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

#ifndef bempp_fmm_transform_hpp
#define bempp_fmm_transform_hpp

#include <vector>

#include "../common/armadillo_fwd.hpp"
#include "../fiber/scalar_traits.hpp"

namespace Bempp
{

// abstract for now, will use Chebyshev as default in future versions
template <typename ValueType>
class FmmTransform
{
public:
	typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

	FmmTransform(unsigned int quadraturePointCount)
	 : 	m_quadraturePoints(3, quadraturePointCount),
		m_quadratureWeights(quadraturePointCount)
	{
	}

	const arma::Col<CoordinateType>& getWeights() const
	{
		return m_quadratureWeights;
	}
	unsigned int quadraturePointCount() const
	{
		return m_quadraturePoints.n_cols;
	}
	arma::Col<CoordinateType> getQuadraturePoint(unsigned int index) const
	{
		return m_quadraturePoints.col(index);
	}

	// multipole to multipole (M2M) translation matrix
	virtual arma::Mat<ValueType> M2M(
		const arma::Col<CoordinateType>& x1, 
		const arma::Col<CoordinateType>& x2) const = 0;

	// multipole to local (M2L) translation matrix
	virtual arma::Mat<ValueType> M2L(
		const arma::Col<CoordinateType>& x1, 
		const arma::Col<CoordinateType>& x2) const = 0;

	// local to local (L2L) translation matrix
	virtual arma::Mat<ValueType> L2L(
		const arma::Col<CoordinateType>& x1, 
		const arma::Col<CoordinateType>& x2) const = 0;

	// multipole expansion coefficients (MEC)
	virtual void evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const = 0;

	// locl expansion coefficients (LEC)
	virtual void evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const = 0;

protected:
	arma::Mat<CoordinateType> m_quadraturePoints;
	arma::Col<CoordinateType> m_quadratureWeights;
};

// all operations are diagonal in the case of plane wave expansion
template <typename ValueType>
class FmmHighFreq : public FmmTransform<ValueType>
{
public:
	typedef typename FmmTransform<ValueType>::CoordinateType CoordinateType;

	FmmHighFreq(ValueType kappa, unsigned int L)//, unsigned int numQuadPoints)
	 : m_kappa(kappa), m_L(L), FmmTransform<ValueType>(L*(2*L+1))//numQuadPoints)//(770)//86)//
	{
		generateGaussPoints();
		// finally all levels in the octree will have their own set of Gauss points
		// will need to interpolate functions here, between each level
		// also cache T matrix
	}

	// multipole to multipole (M2M) translation matrix
	virtual arma::Mat<ValueType> M2M(
		const arma::Col<CoordinateType>& x1, 
		const arma::Col<CoordinateType>& x2) const;

	// multipole to local (M2L) translation matrix
	virtual arma::Mat<ValueType> M2L(
		const arma::Col<CoordinateType>& x1, 
		const arma::Col<CoordinateType>& x2) const;

	// local to local (L2L) translation matrix
	virtual arma::Mat<ValueType> L2L(
		const arma::Col<CoordinateType>& x1, 
		const arma::Col<CoordinateType>& x2) const;

	virtual void generateGaussPoints();

	ValueType kappa() const {return m_kappa;}
private:
	ValueType m_kappa;
	unsigned int m_L;
};

} // namespace Bempp

#endif
