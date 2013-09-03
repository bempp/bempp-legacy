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

#include "fmm_cache.hpp"
#include "fmm_transform.hpp"
#include "octree.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>		// std::cout

namespace Bempp
{

template <typename ValueType>
FmmCache<ValueType>::FmmCache(
		const FmmTransform<ValueType>& fmmTransform,
		unsigned int levels,
		const arma::Col<CoordinateType> &lowerBound,
		const arma::Col<CoordinateType> &upperBound)
	: m_topLevel(2)
{
	arma::Col<CoordinateType> origin(3);
	origin.fill(0.);

	m_cacheM2L.resize(levels-m_topLevel+1);

	for (unsigned int level = m_topLevel; level<=levels; level++) {

		// there are 7^3-3^3=316 unique translation matrices for translation
		// invariant operators
		m_cacheM2L[level-m_topLevel].resize(316);

		unsigned int boxesPerSide = getNodesPerSide(level);
		arma::Col<CoordinateType> boxSize;
		boxSize = (upperBound - lowerBound)/boxesPerSide;

		arma::Col<CoordinateType> centre(3);

		unsigned int index = 0;
		for (int indx=-3; indx<=3; indx++) {
			centre[0] = indx*boxSize[0];
			for (int indy=-3; indy<=3; indy++) {
				centre[1] = indy*boxSize[1];
				for (int indz=-3; indz<=3; indz++) {
					centre[2] = indz*boxSize[2];

					if (abs(indx) > 1 || abs(indy) > 1 || abs(indz) > 1) {

						arma::Mat<ValueType> m2l = fmmTransform.M2L(centre, origin);

						m_cacheM2L[level-m_topLevel][index++] = m2l;

					} // if not a nearest neighbour
				} // for each x offset
			} // for each x offset
		} // for each x offset
	} // for each level

	// M2M & L2L cache
	m_cacheM2M.resize(levels-m_topLevel);
	m_cacheL2L.resize(levels-m_topLevel);
	for (unsigned int level = m_topLevel; level<levels; level++) { //2->levels-1
		m_cacheM2M[level-m_topLevel].resize(8);
		m_cacheL2L[level-m_topLevel].resize(8);

		unsigned int boxesPerSide = getNodesPerSide(level);
		arma::Col<CoordinateType> boxSize;
		boxSize = (upperBound - lowerBound)/boxesPerSide;

		arma::Col<CoordinateType> Rchild(3); // child pos

		for (unsigned long child = 0;	child < 8; child++) {

			unsigned long ind[3];
			deMorton(&ind[0], &ind[1], &ind[2], child);

			Rchild(0) = ind[0]; Rchild(1) = ind[1]; Rchild(2) = ind[2];
			Rchild = (Rchild-0.5) % boxSize/2;

			arma::Mat<ValueType> m2m = fmmTransform.M2M(Rchild, origin);

			m_cacheM2M[level-m_topLevel][child] = m2m;

			arma::Mat<ValueType> l2l = fmmTransform.L2L(origin, Rchild);

			m_cacheL2L[level-m_topLevel][child] = l2l;
		} // for each child
	} // for each level

} // FmmCache

template <typename ValueType>
arma::Mat<ValueType> 
FmmCache<ValueType>::M2M(unsigned int level, unsigned int item) const
{
	return m_cacheM2M[level-m_topLevel][item];
}

template <typename ValueType>
arma::Mat<ValueType> 
FmmCache<ValueType>::M2L(unsigned int level, unsigned int item) const
{
	return m_cacheM2L[level-m_topLevel][item];
}

template <typename ValueType>
arma::Mat<ValueType> 
FmmCache<ValueType>::L2L(unsigned int level, unsigned int item) const
{
	return m_cacheL2L[level-m_topLevel][item];
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmCache);

} // namespace Bempp
