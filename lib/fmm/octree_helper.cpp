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

#include "octree_helper.hpp"
#include "fmm_transform.hpp"
#include "transmit_receive.hpp"
#include "octree.hpp"
#include "octree_node.hpp"

#include "../common/common.hpp"
#include "../common/types.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include "../fiber/types.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/conjugate.hpp"

#include "../space/space.hpp"
#include "../assembly/local_dof_lists_cache.hpp"
#include "../assembly/surface_normal_independent_function.hpp"
#include "../assembly/surface_normal_dependent_function.hpp"

#include "transmit_receive.hpp"
#include "fmm_local_assembler.hpp"

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#include <iostream>		// std::cout
#include <vector>		// std::vector
#include <math.h>		// floor

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
OctreeNearHelper<BasisFunctionType, ResultType>::OctreeNearHelper(
			const shared_ptr<Octree<ResultType> > octree,
			const Space<BasisFunctionType>& testSpace,
			const Space<BasisFunctionType>& trialSpace,
			const std::vector<LocalAssembler*>& assemblers,
               const std::vector<ResultType>& denseTermsMultipliers,
			const AssemblyOptions& options,
			const std::vector<unsigned int> &test_p2o,
			const std::vector<unsigned int> &trial_p2o,
			bool indexWithGlobalDofs)
	:	m_octree(octree), m_testSpace(testSpace), m_trialSpace(trialSpace),
		m_assemblers(assemblers), m_denseTermsMultipliers(denseTermsMultipliers),
		m_options(options)
{
	// requires p2o the vector, not an IndexPermutation object. Provide access functions to cache for nodes
	m_testDofListsCache = boost::make_shared<LocalDofListsCache<BasisFunctionType> >
		(m_testSpace, test_p2o, indexWithGlobalDofs);
	m_trialDofListsCache = boost::make_shared<LocalDofListsCache<BasisFunctionType> >
		(m_trialSpace, trial_p2o, indexWithGlobalDofs);
}



// calculate the near field integrals acting between pairs of octree nodes
// modified from WeakFormAcaAssemblyHelper
// TODO: Add back in Sparse terms
// check the matrix is being returned in an efficient way
// check elements are being evaluated efficiently and not multiple times
template <typename BasisFunctionType, typename ResultType>
arma::Mat<ResultType>
OctreeNearHelper<BasisFunctionType, ResultType>::evaluateNearFieldIntegrals(
		unsigned int dofStartTest, unsigned int dofCountTest,
		unsigned int dofStartTrial, unsigned int dofCountTrial) const
{
	// Convert matrix indices into DOF indices
	shared_ptr<const LocalDofLists<BasisFunctionType> > 
		testDofLists  = m_testDofListsCache->get(
		dofStartTest, dofCountTest);
	shared_ptr<const LocalDofLists<BasisFunctionType> > 
		trialDofLists = m_trialDofListsCache->get(
		dofStartTrial, dofCountTrial);
	// Necessary elements
	const std::vector<int>& testElementIndices  = testDofLists->elementIndices;
	const std::vector<int>& trialElementIndices = trialDofLists->elementIndices;
	// Necessary local dof indices in each element
	const std::vector<std::vector<LocalDofIndex> >& testLocalDofs =
		testDofLists->localDofIndices;
	const std::vector<std::vector<LocalDofIndex> >& trialLocalDofs =
		trialDofLists->localDofIndices;
    // Weights of local dofs in each element
    const std::vector<std::vector<BasisFunctionType> >& testLocalDofWeights =
            testDofLists->localDofWeights;
    const std::vector<std::vector<BasisFunctionType> >& trialLocalDofWeights =
            trialDofLists->localDofWeights;
    for (size_t i = 0; i < testLocalDofWeights.size(); ++i)
        for (size_t j = 0; j < testLocalDofWeights[i].size(); ++j)
            assert(std::abs(testLocalDofWeights[i][j]) > 0.);
    for (size_t i = 0; i < trialLocalDofWeights.size(); ++i)
        for (size_t j = 0; j < trialLocalDofWeights[i].size(); ++j)
            assert(std::abs(trialLocalDofWeights[i][j]) > 0.);
	// Row and column indices in the matrix to be calculated and stored
	const std::vector<std::vector<int> >& blockRows =
		testDofLists->arrayIndices;
	const std::vector<std::vector<int> >& blockCols =
		trialDofLists->arrayIndices;

	arma::Mat<ResultType> result(dofCountTest, dofCountTrial);

	result.fill(0.);

	// First, evaluate the contributions of the dense terms
	// The whole block or its submatrix needed. This means that we are
	// likely to need all or almost all local DOFs from most elements.
	// Evaluate the full local weak form for each pair of test and trial
	// elements and then select the entries that we need.

	Fiber::_2dArray<arma::Mat<ResultType> > localResult;
	for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm) {
		m_assemblers[nTerm]->evaluateLocalWeakForms(
			testElementIndices, trialElementIndices, localResult);//, minDist);
		for (size_t nTrialElem = 0;
			nTrialElem < trialElementIndices.size();
			++nTrialElem)

			for (size_t nTrialDof = 0;
				nTrialDof < trialLocalDofs[nTrialElem].size();
				++nTrialDof)
				for (size_t nTestElem = 0;
					nTestElem < testElementIndices.size();
					++nTestElem)
					for (size_t nTestDof = 0;
						nTestDof < testLocalDofs[nTestElem].size();
						++nTestDof)
						result(blockRows[nTestElem][nTestDof],
							blockCols[nTrialElem][nTrialDof]) +=
							m_denseTermsMultipliers[nTerm] *
							Fiber::conjugate(testLocalDofWeights[nTestElem][nTestDof]) *
							trialLocalDofWeights[nTrialElem][nTrialDof] *
							localResult(nTestElem, nTrialElem)
							(testLocalDofs[nTestElem][nTestDof],
							trialLocalDofs[nTrialElem][nTrialDof]);
	} // for each term

	return result;

} // Octree::evaMatluateNearFieldIntegrals

template <typename BasisFunctionType, typename ResultType>
void OctreeNearHelper<BasisFunctionType, ResultType>::evaluateNearField(
	const shared_ptr<Octree<ResultType> > &octree,
	unsigned int n) const
{
	//const unsigned int nLeaves = getNodesPerLevel(octree->levels());

	//for (unsigned int n=0; n<nLeaves; n++) {
		OctreeNode<ResultType> &node = 
			octree->getNode(n, octree->levels());

		if (node.testDofCount()==0) {
			return; //continue;
		}
		const std::vector<unsigned long>& neigbourList = node.neigbourList();

		// test near field interactions by including all other leaves in neighbour list
		//std::vector<unsigned long> neigbourList(nLeaves-1);
		//for (unsigned int k=1; k<nLeaves; k++) {
		//	neigbourList[k-1] = (n+k)%nLeaves;
		//}
		//node.m_neigbourList = neigbourList;

		// size of the neigbourhood plus the current node
		std::vector<arma::Mat<ResultType> > nearFieldMats;
		nearFieldMats.resize(neigbourList.size()+1);

		const unsigned int testDofStart = node.testDofStart();
		const unsigned int testDofCount = node.testDofCount();

		// first entry is the self-interaction
		if (node.trialDofCount()) {
			unsigned int trialDofStart = node.trialDofStart();
			unsigned int trialDofCount = node.trialDofCount();

			nearFieldMats[0] = evaluateNearFieldIntegrals(
				testDofStart, testDofCount,
				trialDofStart, trialDofCount);
		}

		// repeat for the neighbours: trial functions are fixed in the current node
		// test functions are in the neigbourhood
		for (unsigned long neigh = 0; neigh < neigbourList.size(); neigh++) {
			const OctreeNode<ResultType> &nodeneigh = octree->getNodeConst(
				neigbourList[neigh], octree->levels());
			unsigned int trialDofStart = nodeneigh.trialDofStart();
			unsigned int trialDofCount = nodeneigh.trialDofCount();
	
			nearFieldMats[neigh+1] = evaluateNearFieldIntegrals(
				testDofStart, testDofCount, trialDofStart, trialDofCount);
		}

		node.setNearFieldMats(nearFieldMats);

	//} //for each leaf
}


template <typename BasisFunctionType, typename ResultType>
void OctreeNearHelper<BasisFunctionType, ResultType>::operator()(unsigned int index) const
{
	evaluateNearField(m_octree, index);
}


template <typename BasisFunctionType, typename ResultType>
OctreeFarHelper<BasisFunctionType, ResultType>::OctreeFarHelper(
			const shared_ptr<Octree<ResultType> > octree,
			const Space<BasisFunctionType>& testSpace,
			const Space<BasisFunctionType>& trialSpace,
			const AssemblyOptions& options,
			const std::vector<unsigned int> &test_p2o,
			const std::vector<unsigned int> &trial_p2o,
			bool indexWithGlobalDofs,
			const FmmTransform<ResultType> &fmmTransform)
	:	m_octree(octree), m_testSpace(testSpace), m_trialSpace(trialSpace),
		m_options(options), m_fmmTransform(fmmTransform)
{
	// requires p2o the vector, not an IndexPermutation object. Provide access functions to cache for nodes
	m_testDofListsCache = boost::make_shared<LocalDofListsCache<BasisFunctionType> >
		(m_testSpace, test_p2o, indexWithGlobalDofs);
	m_trialDofListsCache = boost::make_shared<LocalDofListsCache<BasisFunctionType> >
		(m_trialSpace, trial_p2o, indexWithGlobalDofs);
}


template <typename BasisFunctionType, typename ResultType>
arma::Mat<ResultType> //Col
OctreeFarHelper<BasisFunctionType, ResultType>::evaluateFarFieldIntegrals(
		FmmLocalAssembler<BasisFunctionType, ResultType>& fmmLocalAssembler,
		const FmmTransform<ResultType>& fmmTransform,
		const arma::Col<CoordinateType> &nodeCentre, 
		const arma::Col<CoordinateType> &nodeSize, 
		unsigned int dofStart, unsigned int dofCount, bool isTest) const
{
	unsigned int multipoleCount = fmmTransform.quadraturePointCount();

	arma::Mat<ResultType> result(multipoleCount, dofCount);
	result.fill(0.);

	// Convert matrix indices into DOF indices
	shared_ptr<const LocalDofLists<BasisFunctionType> > dofLists;
	if (isTest) {
		dofLists = m_testDofListsCache->get(dofStart, dofCount);
	} else {
		dofLists = m_trialDofListsCache->get(dofStart, dofCount);
	}
	// Necessary elements
	const std::vector<int>& elementIndices = dofLists->elementIndices;
	// Necessary local dof indices in each element
	const std::vector<std::vector<LocalDofIndex> >& localDofs =
		dofLists->localDofIndices;
	// Weights of local dofs in each element
	const std::vector<std::vector<BasisFunctionType> >& localDofWeights =
		dofLists->localDofWeights;
	// Row and column indices in the matrix to be calculated and stored
	const std::vector<std::vector<int> >& blockCols =
		dofLists->arrayIndices;

	// First, evaluate the contributions of the dense terms
	// The whole block or its submatrix needed. This means that we are
	// likely to need all or almost all local DOFs from most elements.
	// Evaluate the full local weak form for each pair of test and trial
	// elements and then select the entries that we need.

	for (size_t multipole = 0; multipole < multipoleCount; ++multipole) {

		arma::Col<CoordinateType> khat = fmmTransform.getQuadraturePoint(multipole);

		typedef ResultType UserFunctionType;

		if (isTest) {
			// functor reference used in function, so must keep in memory
			typedef FmmFunctionMultiplyingTest<UserFunctionType> FunctorType;
			FunctorType functor(khat, nodeCentre, nodeSize, fmmTransform);
			SurfaceNormalDependentFunction<FunctorType> function(functor);

			// duplicated code due to scoping
			fmmLocalAssembler.setFunction(&function);

			std::vector<arma::Col<ResultType> > localResult;
			fmmLocalAssembler.evaluateLocalWeakForms(
				elementIndices, localResult);

			for (size_t nElem = 0; nElem < elementIndices.size(); ++nElem)
				for (size_t nDof = 0;
					nDof < localDofs[nElem].size();
					++nDof)
					result(multipole, blockCols[nElem][nDof]) +=
						Fiber::conjugate(localDofWeights[nElem][nDof])*
						localResult[nElem]
						(localDofs[nElem][nDof]);
		} else {
			typedef FmmFunctionMultiplyingTrial<UserFunctionType> FunctorType;
			FunctorType functor(khat, nodeCentre, nodeSize, fmmTransform);
			SurfaceNormalDependentFunction<FunctorType> function(functor);

			// duplicated code due to scoping
			fmmLocalAssembler.setFunction(&function);

			std::vector<arma::Col<ResultType> > localResult;
			fmmLocalAssembler.evaluateLocalWeakForms(
				elementIndices, localResult);

			for (size_t nElem = 0; nElem < elementIndices.size(); ++nElem)
				for (size_t nDof = 0;
					nDof < localDofs[nElem].size();
					++nDof)
					result(multipole, blockCols[nElem][nDof]) +=
						localDofWeights[nElem][nDof]*
						localResult[nElem]
						(localDofs[nElem][nDof]);
		} // if test

	} // for each multipole

	return result;

} // Octree::evaluateFarFieldMultipoleIntegrals


template <typename BasisFunctionType, typename ResultType>
void OctreeFarHelper<BasisFunctionType, ResultType>::operator()(
	const tbb::blocked_range<unsigned int>& range) const
{
	// will need one FmmLocalAssembler per process, since each leaf modifies
	// function, which is a member function
	FmmLocalAssembler<BasisFunctionType, ResultType> 
		fmmTestLocalAssembler(m_testSpace, m_options, true);

	FmmLocalAssembler<BasisFunctionType, ResultType> 
		fmmTrialLocalAssembler(m_trialSpace, m_options, false);

	for( unsigned int n=range.begin(); n!=range.end(); ++n ) {

		OctreeNode<ResultType> &node = 
			m_octree->getNode(n, m_octree->levels());
	
		arma::Col<CoordinateType> nodeCentre, nodeSize;
		m_octree->nodeCentre(n, m_octree->levels(), nodeCentre);
		m_octree->nodeSize(m_octree->levels(), nodeSize);

		// evaulate and store test far-field
		if (node.testDofCount()) {
			const unsigned int testDofStart = node.testDofStart();
			const unsigned int testDofCount = node.testDofCount();
			arma::Mat<ResultType> testFarFieldMat;
			testFarFieldMat = evaluateFarFieldIntegrals(
				fmmTestLocalAssembler, m_fmmTransform, nodeCentre, 
				nodeSize, testDofStart, testDofCount, true);
			testFarFieldMat = testFarFieldMat.st(); // transpose

			node.setTestFarFieldMat(testFarFieldMat);
		}

		// evaluate and store trial far-field
		if (node.trialDofCount()) {
			const unsigned int trialDofStart = node.trialDofStart();
			const unsigned int trialDofCount = node.trialDofCount();
			arma::Mat<ResultType> trialFarFieldMat;
			trialFarFieldMat = evaluateFarFieldIntegrals(
				fmmTrialLocalAssembler, m_fmmTransform, nodeCentre, 
				nodeSize, trialDofStart, trialDofCount, false);

			node.setTrialFarFieldMat(trialFarFieldMat);
		}

	} // for each node
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(OctreeNearHelper);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(OctreeFarHelper);

} // namespace Bempp
