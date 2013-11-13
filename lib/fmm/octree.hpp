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

#ifndef bempp_octree_hpp
#define bempp_octree_hpp

#include <vector>
#include <complex>

#include "../common/armadillo_fwd.hpp"
#include "../common/common.hpp"
#include "../common/types.hpp"
#include "../assembly/transposition_mode.hpp"

#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"

namespace Bempp
{

unsigned long morton(unsigned long x, unsigned long y, unsigned long z);

void deMorton(	unsigned long *indx, unsigned long *indy, unsigned long *indz, 
		unsigned long n);

unsigned long getParent(unsigned long n);
unsigned long getFirstChild(unsigned long n);
unsigned long getLastChild(unsigned long n);
unsigned long getNodesPerSide(unsigned long level);
unsigned long getNodesPerLevel(unsigned long level);


/** \cond FORWARD_DECL */
template <typename ResultType> class OctreeNode;
template <typename ResultType> class FmmCache;
template <typename ValueType> class FmmTransform;
class IndexPermutation;
/** \endcond */


template <typename ResultType>
class Octree
{
public:
	typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

	Octree(unsigned int levels, 
		const FmmTransform<ResultType> &fmmTransform,
		const shared_ptr<FmmCache<ResultType> > &fmmCache,
		const arma::Col<CoordinateType> &lowerBound,
		const arma::Col<CoordinateType> &upperBound);

	const OctreeNode<ResultType> &getNodeConst(
		unsigned long number, unsigned int level) const;

	void assignPoints(
		bool hermitian, const std::vector<Point3D<CoordinateType> > &testDofCenters,
		const std::vector<Point3D<CoordinateType> > &trialDofCenters,
		std::vector<unsigned int> &test_p2o, std::vector<unsigned int> &trial_p2o);

	void upwardsStep(const FmmTransform<ResultType> &fmmTransform);
	void translationStep(const FmmTransform<ResultType> &fmmTransform);
	void downwardsStep(const FmmTransform<ResultType> &fmmTransform);

	// affects the local and multipole coefficients in the the leaves
	void matrixVectorProduct(
		const arma::Col<ResultType>& x_in,
		arma::Col<ResultType>& y_out,
		const TranspositionMode trans);

	unsigned int levels() const {return m_levels;}
	OctreeNode<ResultType> &getNode(unsigned long number, unsigned int level);
	void nodeCentre(unsigned long number, unsigned int level,
		arma::Col<CoordinateType> &centre) const;
	void nodeSize(unsigned int level,
		arma::Col<CoordinateType> &size) const;
	const FmmCache<ResultType>& fmmCache() {return *m_fmmCache;}
private:
	unsigned long getLeafContainingPoint(const Point3D<CoordinateType> &point) const;
	const unsigned int m_levels;
	const unsigned int m_topLevel;
	// for now use a flat structure
	std::vector<std::vector<OctreeNode<ResultType> > > m_OctreeNodes;
	shared_ptr<IndexPermutation> m_test_p2o, m_trial_p2o;
	const FmmTransform<ResultType>& m_fmmTransform;
	arma::Col<CoordinateType> m_lowerBound, m_upperBound;
	shared_ptr<FmmCache<ResultType> > m_fmmCache;
};

} // namespace Bempp

#endif
