// Copyright (C) 2013 by the BEM++ Authors
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

#ifndef bempp_common_tests_for_spaces_hpp
#define bempp_common_tests_for_spaces_hpp

#include "../check_arrays_are_close.hpp"

#include "common/acc.hpp"
#include "common/scalar_traits.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/boundary_operator.hpp"
#include "assembly/context.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/identity_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"

#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/index_set.hpp"
#include "grid/grid.hpp"
#include "grid/grid_view.hpp"
#include "space/space.hpp"

#include <boost/test/unit_test.hpp>

namespace Bempp
{

template <typename BasisFunctionType>
void local2global_matches_global2local(const Space<BasisFunctionType>& space)
{
    std::unique_ptr<GridView> view = space.grid()->leafView();
    const IndexSet& indexSet = view->indexSet();
    std::vector<int> gdofs;
    std::vector<BasisFunctionType> gdofWeights;
    std::vector<std::vector<LocalDof> > ldofs;
    std::vector<std::vector<BasisFunctionType> > ldofWeights;
    std::unique_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        int elementIndex = indexSet.entityIndex(element);
        space.getGlobalDofs(element, gdofs, gdofWeights);
        BOOST_CHECK_EQUAL(gdofs.size(), gdofWeights.size());
        // some ldofs might be inactive (gdof == -1); to avoid getting an error
        // in the subsequent call to global2localDofs, we pretend that these
        // ldofs correspond to gdof #0.
        for (int i = 0; i < gdofs.size(); ++i)
            if (acc(gdofs, i) < 0)
                acc(gdofs, i) = 0;

        space.global2localDofs(gdofs, ldofs, ldofWeights);
        BOOST_CHECK_EQUAL(ldofs.size(), gdofs.size());
        BOOST_CHECK_EQUAL(ldofs.size(), ldofWeights.size());

        for (int i = 0; i < gdofs.size(); ++i) {
            int gdof = acc(gdofs, i);
            if (gdof == 0) // might in fact be -1 originally, so ignore it
                continue;
            int ldofIndex = -1;
            for (int j = 0; j < acc(ldofs, i).size(); ++j)
                if (acc(ldofs[i], j).entityIndex == elementIndex &&
                    acc(ldofs[i], j).dofIndex == i) {
                    ldofIndex = j;
                    break;
                }
            bool found = ldofIndex >= 0;
            BOOST_CHECK(found == (gdof >= 0));
            if (found)
                BOOST_CHECK(acc(gdofWeights, i) == acc(acc(ldofWeights, i), ldofIndex));
        }

        it->next();
    }
}

template <typename BasisFunctionType>
void global2local_matches_local2global(const Space<BasisFunctionType>& space)
{
    std::unique_ptr<GridView> view = space.grid()->leafView();
    const IndexSet& indexSet = view->indexSet();
    std::vector<std::vector<int> > gdofs; // [elementIndex][ldof]
    std::vector<std::vector<BasisFunctionType> > gdofWeights;
    std::vector<std::vector<LocalDof> > ldofs; // [gdof]
    std::vector<std::vector<BasisFunctionType> > ldofWeights;

    // Fill the arrays gdofs and gdofWeights
    std::unique_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    const int elementCount = view->entityCount(0);
    gdofs.resize(elementCount);
    gdofWeights.resize(elementCount);
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        int elementIndex = indexSet.entityIndex(element);
        space.getGlobalDofs(element, acc(gdofs, elementIndex),
                            acc(gdofWeights, elementIndex));
        BOOST_CHECK_EQUAL(acc(gdofs, elementIndex).size(),
                          acc(gdofWeights, elementIndex).size());
        it->next();
    }

    // Fill the arrays ldofs and ldofWeights
    const int gdofCount = space.globalDofCount();
    std::vector<int> gdofIndices(gdofCount);
    for (int i = 0; i < gdofCount; ++i)
        acc(gdofIndices, i) = i;
    space.global2localDofs(gdofIndices, ldofs, ldofWeights);

    // Check consistency
    for (int gdof = 0; gdof < gdofCount; ++gdof) {
        const std::vector<LocalDof>& activeLdofs = acc(ldofs, gdof);
        const std::vector<BasisFunctionType>& activeLdofWeights =
            acc(ldofWeights, gdof);
        BOOST_CHECK_EQUAL(activeLdofs.size(), activeLdofWeights.size());
        for (int ldofIndex = 0; ldofIndex < activeLdofs.size(); ++ldofIndex) {
            BOOST_CHECK_EQUAL(gdof,
                              acc(acc(gdofs,
                                      acc(activeLdofs, ldofIndex).entityIndex),
                                  acc(activeLdofs, ldofIndex).dofIndex));
            BOOST_CHECK_EQUAL(acc(activeLdofWeights, ldofIndex),
                              acc(acc(gdofWeights,
                                      acc(activeLdofs, ldofIndex).entityIndex),
                                  acc(activeLdofs, ldofIndex).dofIndex));
        }
    }
}

template <typename BasisFunctionType>
void complement_is_really_a_complement(
        const shared_ptr<Space<BasisFunctionType> >& space,
        const shared_ptr<Space<BasisFunctionType> >& space1,
        const shared_ptr<Space<BasisFunctionType> >& space2)
{
    typedef BasisFunctionType BFT;
    typedef BasisFunctionType RT;
    typedef typename ScalarTraits<RT>::RealType CT;

    BOOST_CHECK(space->globalDofCount() > 0);
    BOOST_CHECK(space1->globalDofCount() > 0);
    BOOST_CHECK(space2->globalDofCount() > 0);
    BOOST_CHECK_EQUAL(space->globalDofCount(), 
                      space1->globalDofCount() + space2->globalDofCount());

    AccuracyOptions accuracyOptions;
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));
    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    GridFunction<BFT, RT> gf1(
                context, space1,
                arma::ones<arma::Col<RT> >(space1->globalDofCount()));
    GridFunction<BFT, RT> gf2(
                context, space2,
                arma::ones<arma::Col<RT> >(space2->globalDofCount()));

    BoundaryOperator<BFT, RT> s1_to_s =
            identityOperator<BFT, RT>(context, space1, space, space);
    BoundaryOperator<BFT, RT> s2_to_s =
            identityOperator<BFT, RT>(context, space2, space, space);

    GridFunction<BFT, RT> total = s1_to_s * gf1 + s2_to_s * gf2;
    arma::Col<RT> ones(space->globalDofCount());
    ones.fill(1.);
    BOOST_CHECK(check_arrays_are_close<RT>(total.coefficients(), ones,
                                           100. * std::numeric_limits<CT>::epsilon()));
}

} // namespace Bempp

#endif
