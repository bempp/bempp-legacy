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

#ifndef bempp_fmm_cache_hpp
#define bempp_fmm_cache_hpp

#include <vector>

#include "../common/armadillo_fwd.hpp"
#include "../fiber/scalar_traits.hpp"

namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename ResultType> class FmmTransform;
/** \endcond */

template <typename ValueType>
class FmmCache
{
public:
	typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

	FmmCache(
		const FmmTransform<ValueType>& fmmTransform,
		unsigned int levels,
		const arma::Col<CoordinateType> &lowerBound,
		const arma::Col<CoordinateType> &upperBound);

	arma::Mat<ValueType> M2M(unsigned int level, unsigned int item) const;
	arma::Mat<ValueType> M2L(unsigned int level, unsigned int item) const;
	arma::Mat<ValueType> L2L(unsigned int level, unsigned int item) const;

private:
	const unsigned int m_topLevel; // top level in the octee

	std::vector<std::vector<arma::Mat<ValueType> > > m_cacheM2M;
	std::vector<std::vector<arma::Mat<ValueType> > > m_cacheM2L;
	std::vector<std::vector<arma::Mat<ValueType> > > m_cacheL2L;
};

} // namespace Bempp

#endif
