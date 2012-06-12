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

#include "mass_matrix_container_initialiser.hpp"

#include "config_trilinos.hpp"

#include "mass_matrix_container.hpp"
#include "space.hpp"

#include "../assembly/assembly_options.hpp"
#include "../assembly/discrete_dense_linear_operator.hpp"
#include "../assembly/discrete_inverse_sparse_linear_operator.hpp"
#include "../assembly/discrete_sparse_linear_operator.hpp"
#include "../assembly/identity_operator.hpp"
#include "../assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
MassMatrixContainerInitialiser<BasisFunctionType, ResultType>::
~MassMatrixContainerInitialiser()
{
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<MassMatrixContainer<ResultType> >
MassMatrixContainerInitialiser<BasisFunctionType, ResultType>::operator()() const
{
    typedef DiscreteSparseLinearOperator<ResultType> DiscreteSparseLinOp;
#ifdef WITH_TRILINOS
    typedef DiscreteInverseSparseLinearOperator<ResultType> DiscreteInverseSparseLinOp;
#else
    typedef DiscreteDenseLinearOperator<ResultType> DiscreteDenseLinOp;
#endif

    std::auto_ptr<MassMatrixContainer<ResultType> > result(
                new MassMatrixContainer<ResultType>);
    AssemblyOptions assemblyOptions;
#ifdef WITH_TRILINOS
    // This effectively switches to sparse mode
    assemblyOptions.switchToAca(AcaOptions());
#endif
    IdentityOperator<BasisFunctionType, ResultType> id(m_space, m_space);
    StandardLocalAssemblerFactoryForOperatorsOnSurfaces<
            BasisFunctionType, ResultType> factory;
    result->massMatrix = id.assembleDetachedWeakForm(factory, assemblyOptions);

#ifdef WITH_TRILINOS
    DiscreteSparseLinOp& sparseDiscreteId =
            dynamic_cast<DiscreteSparseLinOp&>(*result->massMatrix);
    Epetra_CrsMatrix& epetraMat = sparseDiscreteId.epetraMatrix();

    result->inverseMassMatrix.reset(
                new DiscreteInverseSparseLinOp(epetraMat, true /* symmetric */));
#else
    // TODO: store an LU decomposition instead of the inverse matrix.
    // (the problem will disappear if we make Trilinos a mandatory dependence).
    result->inverseMassMatrix.reset(
                new DiscreteDenseLinOp(arma::inv(result->massMatrix->asMatrix())));
#endif
    return result;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(MassMatrixContainerInitialiser);

} // namespace Bempp
