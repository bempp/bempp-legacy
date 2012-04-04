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

#include "source_term.hpp"

#include "assembly_options.hpp"
#include "discrete_scalar_valued_source_term.hpp"

#include "../fiber/basis.hpp"
#include "../fiber/local_assembler_factory.hpp"
#include "../fiber/local_assembler_for_source_terms.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../space/space.hpp"

namespace Bempp
{

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedSourceTerm<ValueType> >
SourceTerm<ValueType>::assembleWeakForm(
        const Fiber::Function<ValueType>& function,
        const Space<ValueType>& testSpace,
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    if (!testSpace.dofsAssigned())
        throw std::runtime_error("SourceTerm::assembleWeakForm(): "
                                 "degrees of freedom must be assigned "
                                 "before calling assembleWeakForm()");

    // Prepare local assembler

    std::auto_ptr<GridView> view = testSpace.grid().leafView();
    const int elementCount = view->entityCount(0);

    // Gather geometric data
    Fiber::RawGridGeometry<ValueType> rawGeometry;
    view->getRawElementData(
                rawGeometry.vertices(), rawGeometry.elementCornerIndices(),
                rawGeometry.auxData());

    // Make geometry factory
    std::auto_ptr<GeometryFactory> geometryFactory =
            testSpace.grid().elementGeometryFactory();

    // Get pointers to test and trial bases of each element
    std::vector<const Fiber::Basis<ValueType>*> testBases;
    testBases.reserve(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        testBases.push_back(&testSpace.basis(element));
        it->next();
    }

    // Get reference to the test expression
    const Fiber::Expression<ValueType>& testExpression =
            testSpace.shapeFunctionValueExpression();

    // Now create the assembler
    Fiber::OpenClHandler<ValueType,int> openClHandler(options.openClOptions());
    openClHandler.pushGeometry (rawGeometry.vertices(),
				rawGeometry.elementCornerIndices());

    std::auto_ptr<LocalAssembler> assembler =
            factory.make(*geometryFactory, rawGeometry,
                         testBases,
                         testExpression, function,
                         openClHandler);

    return reallyAssembleWeakForm(
                testSpace, *assembler, options);
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedSourceTerm<ValueType> >
SourceTerm<ValueType>::reallyAssembleWeakForm(
        const Space<ValueType>& testSpace,
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    // Get the grid's leaf view so that we can iterate over elements
    std::auto_ptr<GridView> view = testSpace.grid().leafView();
    const int elementCount = view->entityCount(0);

    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > testGlobalDofs(elementCount);

    // Gather global DOF lists
    const Mapper& mapper = view->elementMapper();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        const int elementIndex = mapper.entityIndex(element);
        testSpace.globalDofs(element, testGlobalDofs[elementIndex]);
        it->next();
    }

    // Make a vector of all element indices
    std::vector<int> testIndices(elementCount);
    for (int i = 0; i < elementCount; ++i)
        testIndices[i] = i;

    // Create the weak form's column vector
    arma::Col<ValueType> result(testSpace.globalDofCount());
    result.fill(0.);

    std::vector<arma::Col<ValueType> > localResult;
    // Evaluate local weak forms
    assembler.evaluateLocalWeakForms(testIndices, localResult);

    // Loop over test indices
    for (int testIndex = 0; testIndex < elementCount; ++testIndex)
        // Add the integrals to appropriate entries in the global weak form
        for (int testDof = 0; testDof < testGlobalDofs[testIndex].size(); ++testDof)
            result(testGlobalDofs[testIndex][testDof]) +=
                    localResult[testIndex](testDof);

    // Create and return a discrete source term represented by the vector that
    // has just been calculated
    return std::auto_ptr<DiscreteScalarValuedSourceTerm<ValueType> >(
                new DiscreteScalarValuedSourceTerm<ValueType>(result));
}


#ifdef COMPILE_FOR_FLOAT
template class SourceTerm<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class SourceTerm<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class SourceTerm<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class SourceTerm<std::complex<double> >;
#endif

} // namespace Bempp
