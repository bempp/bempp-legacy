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

    result->inverseMassMatrix.reset(new DiscreteInverseSparseLinOp(epetraMat));
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
