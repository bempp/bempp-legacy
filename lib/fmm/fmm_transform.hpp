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


#include "../common/armadillo_fwd.hpp"
#include "../fiber/scalar_traits.hpp"

namespace Bempp
{

template <typename ValueType>
ValueType getI();

// abstract for now, will use Chebyshev as default in future versions
// might want to cache translation matrices in future -> store in octree somehow
template <typename ValueType>
class FmmTransform
{
public:
	typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

	FmmTransform(unsigned int P) : m_P(P), m_w(P)
	{
		m_s = new CoordinateType[3*m_P]; // {x,y,z}
//		m_w = new CoordinateType[m_P];
	}
	virtual ~FmmTransform()
	{
		delete[] m_s;
//		delete[] m_w;
	}
	CoordinateType w(unsigned int n) const { return m_w[n]; }
	const arma::Col<CoordinateType>& getWeights() const { return m_w; }

	unsigned int P() const { return m_P; }
	arma::Col<CoordinateType> s(unsigned int ind) const
	{
		arma::Col<CoordinateType> tmp(3);
		tmp(0) = m_s[3*ind+0];
		tmp(1) = m_s[3*ind+1];
		tmp(2) = m_s[3*ind+2];
		return tmp;
	}

	// multipole to multipole (M2M) translation matrix
	virtual arma::Col<ValueType> M2M(CoordinateType x1[3], CoordinateType x2[3]) const = 0;
	// multipole to local (M2L) translation matrix
	virtual arma::Col<ValueType> M2L(CoordinateType x1[3], CoordinateType x2[3]) const = 0;
	// local to local (L2L) translation matrix
	virtual arma::Col<ValueType> L2L(CoordinateType x1[3], CoordinateType x2[3]) const = 0;

	virtual void evaluateTrial( // multipole expansion coefficients (MEC)
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& centre,
			arma::Col<ValueType>& result) const = 0;

	virtual void evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& centre,
			arma::Col<ValueType>& result) const = 0;

protected:
	unsigned int m_P; // number of Gauss points
	CoordinateType *m_s; // Gauss point position
	arma::Col<CoordinateType> m_w; // Gause point weight
};

// all operations are diagonal in the case of plane wave expansion
// TODO: cache some T matrix in particular
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
	virtual arma::Col<ValueType> M2M(CoordinateType x1[3], CoordinateType x2[3]) const;
	// multipole to local (M2L) translation matrix
	virtual arma::Col<ValueType> M2L(CoordinateType x1[3], CoordinateType x2[3]) const;
	// local to local (L2L) translation matrix
	virtual arma::Col<ValueType> L2L(CoordinateType x1[3], CoordinateType x2[3]) const;

	virtual void generateGaussPoints();
	
	ValueType kappa() const {return m_kappa;}
private:
	// must be purely real/imag for now, since boost special functions
	// do not support complex arguments
	ValueType m_kappa;
	unsigned int m_L;
};

} // namespace Bempp

#endif
