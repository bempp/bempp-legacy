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

#include "bempp/common/config_ahmed.hpp"
#include "bempp/common/config_trilinos.hpp"

#include "fmm_global_assembler.hpp"
#include "fmm_transform.hpp"
#include "fmm_cache.hpp"
#include "fmm_far_field_helper.hpp"
#include "fmm_near_field_helper.hpp"
#include "octree.hpp"
#include "octree_node.hpp"
#include "discrete_fmm_boundary_operator.hpp"

#include "../assembly/context.hpp"
//#include "../assembly/evaluation_options.hpp"
#include "../assembly/assembly_options.hpp"
#include "../assembly/index_permutation.hpp"
#include "../assembly/discrete_boundary_operator_composition.hpp"
#include "../assembly/discrete_sparse_boundary_operator.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/auto_timer.hpp"
#include "../common/boost_shared_array_fwd.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../common/chunk_statistics.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/serial_blas_region.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../space/space.hpp"
#include "../grid/grid.hpp"

#include <stdexcept>
#include <iostream>

#include <boost/type_traits/is_complex.hpp>

#include <tbb/atomic.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_queue.h>

#include <complex>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>


namespace Bempp
{

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
FmmGlobalAssembler<BasisFunctionType, KernelType, ResultType>::assembleDetachedWeakForm(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        const std::vector<LocalAssembler*>& localAssemblers,
        const std::vector<const DiscreteBndOp*>& sparseTermsToAdd,
        const std::vector<ResultType>& denseTermsMultipliers,
        const std::vector<ResultType>& sparseTermsMultipliers,
        const Context<BasisFunctionType, ResultType>& context,
        bool hermitian,
        const FmmTransform<ResultType>& fmmTransform,
        const CollectionOfKernels& kernels)
{
	const AssemblyOptions& options = context.assemblyOptions();

	const bool verbosityAtLeastDefault =
            (options.verbosityLevel() >= VerbosityLevel::DEFAULT);
	const bool verbosityAtLeastHigh =
            (options.verbosityLevel() >= VerbosityLevel::HIGH);

	const FmmOptions& fmmOptions = options.fmmOptions();

	const bool indexWithGlobalDofs = true;//fmmOptions.globalAssemblyBeforeCompression;

	const size_t testDofCount = indexWithGlobalDofs ?
		testSpace.globalDofCount() : testSpace.flatLocalDofCount();
	const size_t trialDofCount = indexWithGlobalDofs ?
		trialSpace.globalDofCount() : trialSpace.flatLocalDofCount();

	if (hermitian && testDofCount != trialDofCount)
        throw std::invalid_argument("FmmGlobalAssembler::assembleDetachedWeakForm(): "
                                    "you cannot generate a Hermitian weak form "
                                    "using test and trial spaces with different "
                                    "numbers of DOFs");

	// get the bounding box of the test and trial spaces, in order to set 
	// octree size. Octree should encotestDofCountmpass both spaces
	// N.B. will need a bigger octree for evaulated solution on a surface later
	arma::Col<double> lowerBoundTest, upperBoundTest;
	testSpace.grid()->getBoundingBox(lowerBoundTest, upperBoundTest);

	arma::Col<double> lowerBoundTrial, upperBoundTrial;
	trialSpace.grid()->getBoundingBox(lowerBoundTrial, upperBoundTrial);

	// find the min of each row
	arma::Col<double> lowerBound, upperBound;
	lowerBound = arma::min(arma::join_rows(lowerBoundTest, lowerBoundTrial), 1);
	upperBound = arma::max(arma::join_rows(upperBoundTest, upperBoundTrial), 1);

	if (verbosityAtLeastDefault) {
		std::cout << "lower bound = (" << lowerBound[0] << ", ";
		std::cout << lowerBound[1] << ", " << lowerBound[2] << ')' << std::endl;
		std::cout << "upper bound = (" << upperBound[0] << ", ";
		std::cout << upperBound[1] << ", " << upperBound[2] << ')' << std::endl;

		std::cout << "Caching M2M, M2L and L2L operators" << std::endl;
	}
	shared_ptr<FmmCache<ResultType> > fmmCache = 
		boost::make_shared<FmmCache<ResultType> > (fmmTransform, fmmOptions.levels);
	fmmCache->initCache(
			arma::conv_to<arma::Col<CoordinateType> >::from(lowerBound),
			arma::conv_to<arma::Col<CoordinateType> >::from(upperBound));//,
//			kernels);

	// Note that in future the octree will need to store dof's for test
	// and trial spaces individually, if the two differ in order
	unsigned int nLevels = fmmOptions.levels;
	shared_ptr<Octree<ResultType> > octree = 
		boost::make_shared<Octree<ResultType> >(nLevels, fmmTransform, fmmCache, 
			arma::conv_to<arma::Col<CoordinateType> >::from(lowerBound),
			arma::conv_to<arma::Col<CoordinateType> >::from(upperBound),
			verbosityAtLeastDefault, verbosityAtLeastHigh);

	std::vector<Point3D<CoordinateType> > testDofCenters, trialDofCenters;
	if (indexWithGlobalDofs) {
		testSpace.getGlobalDofPositions(testDofCenters);
		trialSpace.getGlobalDofPositions(trialDofCenters);
	} else {
		testSpace.getFlatLocalDofPositions(testDofCenters);
		trialSpace.getFlatLocalDofPositions(trialDofCenters);
	}
	std::vector<unsigned int> trial_p2o, test_p2o;
	octree->assignPoints(hermitian, testDofCenters, trialDofCenters,
		test_p2o, trial_p2o);

	if (verbosityAtLeastDefault) {
		std::cout << "Caching near field interactions" << std::endl;
	}

	FmmNearFieldHelper<BasisFunctionType, ResultType> fmmNearFieldHelper(
		octree, testSpace, trialSpace, localAssemblers, denseTermsMultipliers, 
		options, test_p2o, trial_p2o, indexWithGlobalDofs);

	unsigned int nLeaves = getNodesPerLevel(octree->levels());
	tbb::parallel_for<unsigned int>(0, nLeaves, fmmNearFieldHelper);
	//octreeHelper.evaluateNearField(octree);

	//std::cout << "Caching trial far-field interactions" << std::endl;
	//octreeHelper.evaluateTrialFarField(octree, fmmTransform);

	//std::cout << "Caching test far-field interactions" << std::endl;
	//octreeHelper.evaluateTestFarField(octree, fmmTransform);

	if (verbosityAtLeastDefault) {
		std::cout << "Caching test and trial far-field interactions" << std::endl;
	}

	FmmFarFieldHelper<BasisFunctionType, ResultType> fmmFarFieldHelper(
		octree, testSpace, trialSpace, options, test_p2o, trial_p2o, 
		indexWithGlobalDofs, fmmTransform);

	tbb::parallel_for(tbb::blocked_range<unsigned int>(0, nLeaves, 100), 
		fmmFarFieldHelper);

	unsigned int symmetry = NO_SYMMETRY;
	if (hermitian) {
		symmetry |= HERMITIAN;
		if (boost::is_complex<ResultType>())
			symmetry |= SYMMETRIC;
	}

	// this is a problem, need to hide BasisFunctionType argument somehow
	typedef DiscreteFmmBoundaryOperator<ResultType> DiscreteFmmLinOp;
	std::auto_ptr<DiscreteFmmLinOp> fmmOp(
                new DiscreteFmmLinOp(testDofCount, trialDofCount,
                                     octree,
                                     Symmetry(symmetry) ));

	std::auto_ptr<DiscreteBndOp> result;
	result = fmmOp;

	return result;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
FmmGlobalAssembler<BasisFunctionType, KernelType, ResultType>::assembleDetachedWeakForm(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        LocalAssembler& localAssembler,
        const Context<BasisFunctionType, ResultType>& context,
        bool hermitian,
        const FmmTransform<ResultType>& fmmTransform,
        const CollectionOfKernels& kernels)
{
    std::vector<LocalAssembler*> localAssemblers(1, &localAssembler);
    std::vector<const DiscreteBndOp*> sparseTermsToAdd;
    std::vector<ResultType> denseTermsMultipliers(1, 1.0);
    std::vector<ResultType> sparseTermsMultipliers;

    return assembleDetachedWeakForm(testSpace, trialSpace, localAssemblers,
                            sparseTermsToAdd,
                            denseTermsMultipliers,
                            sparseTermsMultipliers,
                            context, hermitian, fmmTransform, kernels);
}

//FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmGlobalAssembler);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(FmmGlobalAssembler);

} // namespace Bempp
