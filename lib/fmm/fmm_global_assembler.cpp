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
#include "../common/types.hpp" // defines LocalDof
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/serial_blas_region.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../space/space.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"

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

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
FmmGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        const std::vector<LocalAssembler*>& localAssemblers,
        const std::vector<const DiscreteBndOp*>& sparseTermsToAdd,
        const std::vector<ResultType>& denseTermsMultipliers,
        const std::vector<ResultType>& sparseTermsMultipliers,
        const Context<BasisFunctionType, ResultType>& context,
        bool hermitian,
        const FmmTransform<ResultType>& fmmTransform)
{
    const AssemblyOptions& options = context.assemblyOptions();

    const bool verbosityAtLeastDefault =
            (options.verbosityLevel() >= VerbosityLevel::DEFAULT);
    const bool verbosityAtLeastHigh =
            (options.verbosityLevel() >= VerbosityLevel::HIGH);

    const FmmOptions& fmmOptions = options.fmmOptions();

    // cannot use flatlocaldofs with a triangle based octree and with piecewise 
    // linear elements. For triangles near leaf boarders, the nodes will have 
    // contributions from other attached elements both inside and outside the current
    // leaf. For this reason, a given node should be treated as a far field interaction
    // in one test triangle, but as a far-field interation in a different test triangle.
    const bool indexWithGlobalDofs = true;

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
            arma::conv_to<arma::Col<CoordinateType> >::from(upperBound));

    // Note that in future the octree will need to store dof's for test
    // and trial spaces individually, if the two differ in order
    unsigned int nLevels = fmmOptions.levels;
    shared_ptr<Octree<ResultType> > octree = 
        boost::make_shared<Octree<ResultType> >(nLevels, fmmTransform, fmmCache, 
            arma::conv_to<arma::Col<CoordinateType> >::from(lowerBound),
            arma::conv_to<arma::Col<CoordinateType> >::from(upperBound),
            verbosityAtLeastDefault, verbosityAtLeastHigh);

    // assign each dof a location which is used to determine its leaf in the octree
    // using the barycentre of the triangle for discontinuous spaces
    std::vector<Point3D<CoordinateType> > testDofLocations, trialDofLocations;
    if (trialSpace.isDiscontinuous()) { // default triangle indexed octree

        trialDofLocations.resize(trialDofCount);

        // Form a vector containing barycentres of all triangles in the trial grid,
        // following the approach of grid_segment.cpp

        // must use a const reference to the view to call entityIterator
        const GridView &view = trialSpace.gridView();

        std::vector<Point3D<CoordinateType> > barycentres( view.entityCount(0) );

        std::auto_ptr<EntityIterator<0> > it = view.entityIterator<0>();
        const IndexSet& indexSet = view.indexSet();
        arma::Col<double> barycentre;

        while (!it->finished()) {
            const Entity<0>& entity = it->entity();
            entity.geometry().getCenter(barycentre);

            unsigned int element = indexSet.entityIndex(entity);
            barycentres[element].x = barycentre(0);
            barycentres[element].y = barycentre(1);
            barycentres[element].z = barycentre(2);
            it->next();
        }
    
        // Following local_dof_lists_cache.hpp, find owning triangle of a dof
        typedef int DofIndex;
        std::vector<DofIndex> indices(trialDofCount);
        std::iota(indices.begin(), indices.end(), 0); // fill 0..trialDofCount-1
        std::vector<LocalDof> localDofs;
        trialSpace.flatLocal2localDofs(indices, localDofs);

        for (unsigned int dof = 0; dof < trialDofCount; dof++) {
            unsigned int element = localDofs[dof].entityIndex;
            trialDofLocations[dof] = barycentres[element];
        }

    } else {
        trialSpace.getGlobalDofPositions(trialDofLocations);
    }

    // repeat for the test space
    if (testSpace.isDiscontinuous()) { // default triangle indexed octree

        testDofLocations.resize(testDofCount);

        // Form a vector containing barycentres of all triangles in the test grid,
        // following the approach of grid_segment.cpp

        // must use a const reference to the view to call entityIterator
        const GridView &view = testSpace.gridView();

        std::vector<Point3D<CoordinateType> > barycentres( view.entityCount(0) );

        std::auto_ptr<EntityIterator<0> > it = view.entityIterator<0>();
        const IndexSet& indexSet = view.indexSet();
        arma::Col<double> barycentre;

        while (!it->finished()) {
            const Entity<0>& entity = it->entity();
            entity.geometry().getCenter(barycentre);

            unsigned int element = indexSet.entityIndex(entity);
            barycentres[element].x = barycentre(0);
            barycentres[element].y = barycentre(1);
            barycentres[element].z = barycentre(2);
            it->next();
        }

        // Following local_dof_lists_cache.hpp, find owning triangle of a dof
        typedef int DofIndex;
        std::vector<DofIndex> indices(testDofCount);
        std::iota(indices.begin(), indices.end(), 0); // fill 0..testDofCount-1
        std::vector<LocalDof> localDofs;
        testSpace.flatLocal2localDofs(indices, localDofs);

        for (unsigned int dof = 0; dof < testDofCount; dof++) {
            unsigned int element = localDofs[dof].entityIndex;
            testDofLocations[dof] = barycentres[element];
        }

    } else {
        testSpace.getGlobalDofPositions(testDofLocations);
    }

    std::vector<unsigned int> trial_p2o, test_p2o;
    octree->assignPoints(hermitian, testDofLocations, trialDofLocations,
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

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
FmmGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        LocalAssembler& localAssembler,
        const Context<BasisFunctionType, ResultType>& context,
        bool hermitian,
        const FmmTransform<ResultType>& fmmTransform)
{
    std::vector<LocalAssembler*> localAssemblers(1, &localAssembler);
    std::vector<const DiscreteBndOp*> sparseTermsToAdd;
    std::vector<ResultType> denseTermsMultipliers(1, 1.0);
    std::vector<ResultType> sparseTermsMultipliers;

    return assembleDetachedWeakForm(testSpace, trialSpace, localAssemblers,
                            sparseTermsToAdd,
                            denseTermsMultipliers,
                            sparseTermsMultipliers,
                            context, hermitian, fmmTransform);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmGlobalAssembler);

} // namespace Bempp
