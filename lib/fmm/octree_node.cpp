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

#include "octree_node.hpp"
#include "octree.hpp"

#include "../common/common.hpp"
#include "../common/types.hpp"
#include "../fiber/types.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../common/shared_ptr.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include <iostream>		// std::cout
#include <algorithm>	// std::set_difference, std::sort
#include <vector>		// std::vector
#include <math.h>		// floor
#include <functional>   // std::plus
#include <complex>

namespace Bempp
{


template <typename ResultType>
OctreeNode<ResultType>::OctreeNode(unsigned long number, unsigned int level) 
		:	m_number(number), m_level(level), m_testDofStart(0), m_testDofCount(0),
			m_trialDofStart(0), m_trialDofCount(0)
{
}
//, m_neigbourList(26), m_InteractionList(189) {} this messes up push_back
// note that our neighbour lists neglect the current node

// must be a bit careful, since our neighbour lists are stored explicity
// without empty boxes
template <typename ResultType>
void OctreeNode<ResultType>::makeNeigbourList(const Octree<ResultType> &octree) {

	// can of course still be empty and have neighbours, but it turns
	// out we don't need neighbour lists for empty boxes
	if (testDofCount()==0) {
		return;	// get same interaction lists with or without this line
	}
	// Will probably need the full list of neighbours for reduced scheme
	unsigned long ind[3];
	deMorton(&ind[0], &ind[1], &ind[2], number());

	// find min and max indices along x,y,z spanned by the neighbours
	unsigned long indMin[3], indMax[3];
	for (unsigned int d=0; d<3; d++) {
		indMin[d] = ((ind[d]==0) ? 0 : ind[d]-1); // be careful, unsigned!
		indMax[d] = std::min<unsigned int>( getNodesPerSide(level())-1, ind[d]+1 );
	} // for each dimension, x, y, z

	unsigned long indNeigh[3];
	for (		indNeigh[0]=indMin[0]; indNeigh[0]<=indMax[0]; indNeigh[0]++) {
		for (	indNeigh[1]=indMin[1]; indNeigh[1]<=indMax[1]; indNeigh[1]++) {
			for (indNeigh[2]=indMin[2]; indNeigh[2]<=indMax[2]; indNeigh[2]++) {
				if (indNeigh[0]!=ind[0] || indNeigh[1]!=ind[1] || indNeigh[2]!=ind[2]) {
					unsigned long numberNeigh;
					numberNeigh = morton(indNeigh[0], indNeigh[1], indNeigh[2]);
					if (octree.getNodeConst(numberNeigh, level()).trialDofCount()) {
						m_neigbourList.push_back(numberNeigh);
					} // if the neighbour is not empty
				} // if neighbour is not the current node
			} // for z
		} // for y
	} // for x
	std::sort (m_neigbourList.begin(), m_neigbourList.end());
} // OctreeNode::makeNeigbourList

// genererate the interation list. Need access to the Octree, to get references
// to other nodes. We want to check if these nodes are empty.
// call this function after assigning points to the tree
// only pass in octree where we need it, incase it moves
template <typename ResultType>
void OctreeNode<ResultType>::makeInteractionList(const Octree<ResultType> &octree) {

	if (testDofCount()==0) {
		return;
	}

	// isolated points, might not have any immediate neighbours, but
	// might interact on a larger scale 
	//if (m_neigbourList.size()==0) {
	//	throw new NeighbourListMissing();
	//}

	unsigned int parent = getParent(number());

	// for level 2, assume all boxes besides parent are neighbours
	std::vector<unsigned long> parentNeigbourListLevel2(7);
	for (int k=0; k<7; k++) {
		parentNeigbourListLevel2[k] = (parent+k+1)%8;
	}

	const std::vector<unsigned long> &parentNeigbourList 
		= ((level()==2) ? parentNeigbourListLevel2 : 
		octree.getNodeConst(parent, level()-1).m_neigbourList);

	std::vector<unsigned long> childrenOfParentNeighList;
	for(	std::vector<unsigned long>::const_iterator 
		parentNeighbour = parentNeigbourList.begin(); 
		parentNeighbour != parentNeigbourList.end(); ++parentNeighbour) {

		// loop over the Morton indices of the parent's children
		for (unsigned long child = getFirstChild(*parentNeighbour); 
			child <= getLastChild(*parentNeighbour); child++) {

			if (octree.getNodeConst(child, level()).trialDofCount()) {
				childrenOfParentNeighList.push_back(child);
			}
		}
	} // for each neighbour of the parent

	std::sort (childrenOfParentNeighList.begin(), childrenOfParentNeighList.end());
	std::sort (m_neigbourList.begin(), m_neigbourList.end());

	m_InteractionList.resize(childrenOfParentNeighList.size());
	std::vector<unsigned long>::iterator it;
	it = std::set_difference (childrenOfParentNeighList.begin(), childrenOfParentNeighList.end(),
		m_neigbourList.begin(), m_neigbourList.end(), m_InteractionList.begin());
	m_InteractionList.resize(it-m_InteractionList.begin());
	//std::cout << m_InteractionList.size() << ' '; 
	// the parents with children that do not intersect the neigbour list
	// could be detected here
} // OctreeNode::makeInteractionList

// return the centre of the node in (0,1)^3
/*
template <typename ResultType>
void OctreeNode<ResultType>::centre(CoordinateType centre[3]) const
{
	unsigned long ind[3];
	deMorton(&ind[0], &ind[1], &ind[2], number());
	for (unsigned int d=0; d<3; d++) {
		centre[d] = (ind[d] + 0.5)/getNodesPerSide(level())*2.0-1.0;
	}
} // OctreeNode::centre
*/
template <typename ResultType>
void OctreeNode<ResultType>::setIndex(unsigned long number, unsigned int level) {
	m_number = number;
	m_level = level;
}

// return Morton index
template <typename ResultType>
unsigned long OctreeNode<ResultType>::number() const
{
	return m_number;
}

// return level in the octree
template <typename ResultType>
unsigned int OctreeNode<ResultType>::level() const
{
	return m_level;
}

template <typename ResultType>
ResultType OctreeNode<ResultType>::mcoef(unsigned int n) const
{
	return m_mcoef[n];
}


template <typename ResultType>
void OctreeNode<ResultType>::setMultipoleCoefficients(
		const arma::Col<ResultType> &cvec)
{
	m_mcoef = cvec;
}

template <typename ResultType>
ResultType OctreeNode<ResultType>::lcoef(unsigned int n) const
{
	return m_lcoef[n];
}

template <typename ResultType>
void OctreeNode<ResultType>::setLocalCoefficients(
		const arma::Col<ResultType> &cvec)
{
	m_lcoef = cvec;
}


template <typename ResultType>
void OctreeNode<ResultType>::addLocalCoefficients(
		const arma::Col<ResultType> &cvec)
{
	m_lcoef += cvec;
}

template <typename ResultType>
unsigned int OctreeNode<ResultType>::interactionListSize() const
{
	return m_InteractionList.size();
}

template <typename ResultType>
unsigned int OctreeNode<ResultType>::interactionItem(unsigned int n) const {
	return m_InteractionList[n];
}

template <typename ResultType>
void OctreeNode<ResultType>::setTestDofStart(unsigned int start)
{
	m_testDofStart = start;
}

template <typename ResultType>
void OctreeNode<ResultType>::setTrialDofStart(unsigned int start)
{
	m_trialDofStart = start;
}

template <typename ResultType>
unsigned int OctreeNode<ResultType>::postIncTestDofCount()
{
	return m_testDofCount++;
}

template <typename ResultType>
unsigned int OctreeNode<ResultType>::postIncTrialDofCount()
{
	return m_trialDofCount++;
}





template <typename ResultType>
void OctreeNode<ResultType>::evaluateNearFieldMatrixVectorProduct(
		const Octree<ResultType> &octree,
		const arma::Col<ResultType>& x_in,
		arma::Col<ResultType>& y_in_out)
{
	if (testDofCount()==0) {
		return;
	}

	// same trial functions for this leaf's interactions
	const unsigned int testStart = testDofStart();
	const unsigned int testEnd = testStart + testDofCount() - 1;

	// first entry is the self-interaction
	if (trialDofCount()) {
		unsigned int trialStart = trialDofStart();
		unsigned int trialEnd = trialStart + trialDofCount() - 1;

		const arma::Col<ResultType>& xLocal = x_in.rows(trialStart, trialEnd);
		y_in_out.rows(testStart, testEnd) += m_nearFieldMats[0]*xLocal;
	}

	// repeat for the neighbours: trial functions are fixed in the current node
	// test functions are in the neigbourhood
	for (unsigned long neigh = 0; neigh < m_neigbourList.size(); neigh++) {
		const OctreeNode &node = octree.getNodeConst(m_neigbourList[neigh], m_level);

		//if (node.trialDofCount()==0) {
		//	continue;
		//}

		unsigned int trialStart = node.trialDofStart();
		unsigned int trialEnd = trialStart + node.trialDofCount() - 1;

		const arma::Col<ResultType>& xLocal = x_in.rows(trialStart, trialEnd);
		y_in_out.rows(testStart, testEnd) += m_nearFieldMats[neigh+1]*xLocal;
	}
}




template <typename ResultType>
void OctreeNode<ResultType>::evaluateMultipoleCoefficients(
	const arma::Col<ResultType>& x_in)
{
	if (trialDofCount()==0) {
		return;
	}

	const unsigned int trialStart = trialDofStart();
	const unsigned int trialCount = trialStart + trialDofCount() - 1;

	const arma::Col<ResultType>& xLocal = x_in.rows(trialStart, trialCount);

	// m_trialFarFieldMat(multipole coefficients, dofTrial)
	m_mcoef = m_trialFarFieldMat*xLocal;
}

// local coefficients in each leaf, to far field contributation at each test dof
template <typename ResultType>
void OctreeNode<ResultType>::evaluateFarFieldMatrixVectorProduct(
	const arma::Col<CoordinateType>& weights, arma::Col<ResultType>& y_inout)
{
	if (testDofCount()==0) {
		return;
	}

	const unsigned int testStart = testDofStart();
	const unsigned int testCount = testStart + testDofCount() - 1;

	// m_testFarFieldMat(dofTest, local coefficients)
	// ylocal += [ (diag(lcoef)*m_testFarFieldMat^T)^T ]*weight
	// equilvalent to ylocal += [ m_testFarFieldMat*(weight.*lcoef)
	// part between [] is the local coefficient at each test dof, calc weighted sum
	// special case, simplification possible since L2L operation is diagonal
	y_inout.rows(testStart, testCount) += m_testFarFieldMat*(m_lcoef%weights);
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(OctreeNode);

} // namespace Bempp
