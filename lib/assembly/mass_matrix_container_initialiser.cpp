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

#include "mass_matrix_container.hpp"
#include "../space/space.hpp"

#include "../assembly/assembly_options.hpp"
#include "../assembly/discrete_dense_boundary_operator.hpp"
#include "../assembly/identity_operator.hpp"
#include "../assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "../fiber/explicit_instantiation.hpp"

#ifdef WITH_TRILINOS

#include "../assembly/discrete_inverse_sparse_boundary_operator.hpp"
#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../assembly/discrete_boundary_operator_composition.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_LocalMap.h>
#include <Epetra_SerialComm.h>
#include <EpetraExt_MatrixMatrix.h>
#include <boost/make_shared.hpp>

#endif

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
MassMatrixContainerInitialiser<BasisFunctionType, ResultType>::
MassMatrixContainerInitialiser(const Space<BasisFunctionType>& space,
                               const Space<BasisFunctionType>& dualSpace,
                               bool forceDense) :
    m_space(space), m_dualSpace(dualSpace), m_forceDense(forceDense)
{
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<MassMatrixContainer<ResultType> >
MassMatrixContainerInitialiser<BasisFunctionType, ResultType>::operator()() const
{
#ifdef WITH_TRILINOS
    if (m_forceDense)
        return assembleOperatorsInDenseMode();
    else
        return assembleOperatorsInSparseMode();
#else
    return assembleOperatorsInDenseMode();
#endif
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<MassMatrixContainer<ResultType> >
MassMatrixContainerInitialiser<BasisFunctionType, ResultType>::
assembleOperatorsInDenseMode() const
{
    typedef DiscreteDenseBoundaryOperator<ResultType> DiscreteDenseLinOp;
    std::auto_ptr<MassMatrixContainer<ResultType> > result(
                new MassMatrixContainer<ResultType>);

    IdentityOperator<BasisFunctionType, ResultType> id(m_space, m_space, m_dualSpace);
    AssemblyOptions assemblyOptions;
    StandardLocalAssemblerFactoryForOperatorsOnSurfaces<
            BasisFunctionType, ResultType> factory;
    id.assembleWeakForm(factory, assemblyOptions);
    result->massMatrix = id.weakForm();

    if (result->massMatrix->rowCount() == result->massMatrix->columnCount())
        // TODO: store an LU decomposition instead of the inverse matrix.
        // (the problem will disappear if we make Trilinos a mandatory dependence).
        result->massMatrixPseudoinverse.reset(
                    new DiscreteDenseLinOp(
                        arma::inv(result->massMatrix->asMatrix())));
    else
        // compute and store pseudoinverse
        result->massMatrixPseudoinverse.reset(
                    new DiscreteDenseLinOp(
                        arma::pinv(result->massMatrix->asMatrix())));
    return result;
}

#ifdef WITH_TRILINOS
template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<MassMatrixContainer<ResultType> >
MassMatrixContainerInitialiser<BasisFunctionType, ResultType>::
assembleOperatorsInSparseMode() const
{
    typedef DiscreteBoundaryOperator<ResultType> DiscreteLinOp;
    typedef DiscreteSparseBoundaryOperator<ResultType> DiscreteSparseLinOp;
    typedef DiscreteInverseSparseBoundaryOperator<ResultType> DiscreteInverseSparseLinOp;

    std::auto_ptr<MassMatrixContainer<ResultType> > result(
                new MassMatrixContainer<ResultType>);

    // Construct the mass matrix
    IdentityOperator<BasisFunctionType, ResultType> id(m_space, m_space, m_dualSpace);
    StandardLocalAssemblerFactoryForOperatorsOnSurfaces<
            BasisFunctionType, ResultType> factory;
    AssemblyOptions assemblyOptions;
    // This effectively switches to sparse mode
    assemblyOptions.switchToAca(AcaOptions());
    id.assembleWeakForm(factory, assemblyOptions);
    shared_ptr<const DiscreteLinOp> discreteId = id.weakForm();
    result->massMatrix = discreteId;

    // Construct the mass matrix inverse
    shared_ptr<const DiscreteSparseLinOp> discreteSparseId =
            boost::dynamic_pointer_cast<const DiscreteSparseLinOp>(discreteId);
    shared_ptr<const Epetra_CrsMatrix> idMatrix = discreteSparseId->epetraMatrix();

    const int rowCount = idMatrix->NumGlobalRows();
    const int colCount = idMatrix->NumGlobalCols();
    if (rowCount == colCount) {
        // Square matrix; construct M^{-1}
        bool sameSpace = &m_space == &m_dualSpace;
        result->massMatrixPseudoinverse =
                boost::make_shared<DiscreteInverseSparseLinOp>(
                    idMatrix, sameSpace ? HERMITIAN : NO_SYMMETRY);
    } else {
        // Construct the discrete operator representing M^H
        shared_ptr<DiscreteLinOp> transposeOp =
                boost::make_shared<DiscreteSparseLinOp>(
                    idMatrix, NO_SYMMETRY, CONJUGATE_TRANSPOSE);

        // Construct the discrete operator representing the smaller of
        // (M^H M)^{-1} and (M M^H)^{-1}
        Epetra_SerialComm comm; // To be replaced once we begin to use MPI
        int size = std::min(rowCount, colCount);
        Epetra_LocalMap map(size, 0 /* index_base */, comm);
        shared_ptr<Epetra_CrsMatrix> productMatrix =
                boost::make_shared<Epetra_CrsMatrix>(
                    Copy, map, map,
                    // estimate of the number of nonzero entries per row
                    3 * idMatrix->GlobalMaxNumEntries());

        if (rowCount > colCount) {
            // Tall matrix (overdetermined least-square problem);
            // construct (M^H M)^{-1} M^H
            EpetraExt::MatrixMatrix::Multiply(*idMatrix, true /* transpose */,
                                              *idMatrix, false /* no transpose */,
                                              *productMatrix);
            shared_ptr<DiscreteLinOp> productInverseOp =
                    boost::make_shared<DiscreteInverseSparseLinOp>(
                        productMatrix, HERMITIAN);

            result->massMatrixPseudoinverse =
                    boost::make_shared<DiscreteBoundaryOperatorComposition<ResultType> >(
                        productInverseOp, transposeOp);
        } else {
            // Wide matrix (underdetermined least-square problem);
            // construct M^H (M M^H)^{-1}
            EpetraExt::MatrixMatrix::Multiply(*idMatrix, false /* no transpose */,
                                              *idMatrix, true /* transpose */,
                                              *productMatrix);
            shared_ptr<DiscreteLinOp> productInverseOp =
                    boost::make_shared<DiscreteInverseSparseLinOp>(
                        productMatrix, HERMITIAN);
            result->massMatrixPseudoinverse =
                    boost::make_shared<DiscreteBoundaryOperatorComposition<ResultType> >(
                        transposeOp, productInverseOp);
        }
    }

    return result;
}
#endif

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(MassMatrixContainerInitialiser);

} // namespace Bempp
