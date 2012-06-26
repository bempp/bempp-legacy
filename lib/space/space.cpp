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

#include "space.hpp"

//#include "mass_matrix_container.hpp"
//#include "mass_matrix_container_initialiser.hpp"

//#include "../assembly/discrete_boundary_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"

//#include <boost/type_traits/is_same.hpp>
//#include <boost/static_assert.hpp>

namespace Bempp
{

template <typename BasisFunctionType>
Space<BasisFunctionType>::Space(Grid& grid) :
    m_grid(grid)
//  ,
//    m_bftMassMatrixContainer(
//        MassMatrixContainerInitialiser<BasisFunctionType, BasisFunctionType>(*this)),
//    m_ctMassMatrixContainer(
//        MassMatrixContainerInitialiser<BasisFunctionType, ComplexType>(*this))
{
}

template <typename BasisFunctionType>
Space<BasisFunctionType>::~Space()
{
}

//template <typename BasisFunctionType>
//template <typename ResultType>
//typename boost::enable_if_c<
//boost::is_same<ResultType, BasisFunctionType>::value ||
//boost::is_same<ResultType,
//typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>::value
//>::type
//Space<BasisFunctionType>::applyMassMatrix(
//        const arma::Col<ResultType>& argument,
//        arma::Col<ResultType>& result) const
//{
//    // The compile-time type test ensures that reinterpret_cast is a no-op
//    result.set_size(argument.n_rows);
//    if (boost::is_same<ResultType, BasisFunctionType>::value)
//        applyMassMatrixBasisFunctionType(
//                    reinterpret_cast<const arma::Col<BasisFunctionType>&>(argument),
//                    reinterpret_cast<arma::Col<BasisFunctionType>&>(result));
//    else
//        applyMassMatrixComplexType(
//                    reinterpret_cast<const arma::Col<ComplexType>&>(argument),
//                    reinterpret_cast<arma::Col<ComplexType>&>(result));
//}

//template <typename BasisFunctionType>
//void Space<BasisFunctionType>::applyMassMatrixBasisFunctionType(
//        const arma::Col<BasisFunctionType>& argument,
//        arma::Col<BasisFunctionType>& result) const
//{
//    m_bftMassMatrixContainer.get().massMatrix->apply(
//                NO_TRANSPOSE, argument, result, 1., 0.);
//}

//template <typename BasisFunctionType>
//void Space<BasisFunctionType>::applyMassMatrixComplexType(
//        const arma::Col<ComplexType>& argument,
//        arma::Col<ComplexType>& result) const
//{
//    m_ctMassMatrixContainer.get().massMatrix->apply(
//                NO_TRANSPOSE, argument, result, 1., 0.);
//}

//template <typename BasisFunctionType>
//template <typename ResultType>
//typename boost::enable_if_c<
//boost::is_same<ResultType, BasisFunctionType>::value ||
//boost::is_same<ResultType,
//typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>::value
//>::type
//Space<BasisFunctionType>::applyInverseMassMatrix(
//        const arma::Col<ResultType>& argument,
//        arma::Col<ResultType>& result) const
//{
//    // The compile-time type test ensures that reinterpret_cast is a no-op
//    result.set_size(argument.n_rows);
//    if (boost::is_same<ResultType, BasisFunctionType>::value)
//        applyInverseMassMatrixBasisFunctionType(
//                    reinterpret_cast<const arma::Col<BasisFunctionType>&>(argument),
//                    reinterpret_cast<arma::Col<BasisFunctionType>&>(result));
//    else
//        applyInverseMassMatrixComplexType(
//                    reinterpret_cast<const arma::Col<ComplexType>&>(argument),
//                    reinterpret_cast<arma::Col<ComplexType>&>(result));
//}

//template <typename BasisFunctionType>
//void Space<BasisFunctionType>::applyInverseMassMatrixBasisFunctionType(
//        const arma::Col<BasisFunctionType>& argument,
//        arma::Col<BasisFunctionType>& result) const
//{
//    m_bftMassMatrixContainer.get().massMatrixPseudoinverse->apply(
//                NO_TRANSPOSE, argument, result, 1., 0.);
//}

//template <typename BasisFunctionType>
//void Space<BasisFunctionType>::applyInverseMassMatrixComplexType(
//        const arma::Col<ComplexType>& argument,
//        arma::Col<ComplexType>& result) const
//{
//    m_ctMassMatrixContainer.get().massMatrixPseudoinverse->apply(
//                NO_TRANSPOSE, argument, result, 1., 0.);
//}

/** \brief Get pointers to Basis objects corresponding to all elements. */
template <typename BasisFunctionType>
void getAllBases(const Space<BasisFunctionType>& space,
        std::vector<const Fiber::Basis<BasisFunctionType>*>& bases)
{
    std::auto_ptr<GridView> view = space.grid().leafView();
    const Mapper& mapper = view->elementMapper();
    const int elementCount = view->entityCount(0);

    bases.resize(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& e = it->entity();
        bases[mapper.entityIndex(e)] = &space.basis(e);
        it->next();
    }
}

#define INSTANTIATE(TYPE) \
    template \
    void getAllBases(const Space<TYPE>& space, \
                std::vector<const Fiber::Basis<TYPE>*>& bases)

#ifdef ENABLE_SINGLE_PRECISION
INSTANTIATE(float);
#  ifdef ENABLE_COMPLEX_BASIS_FUNCTIONS
INSTANTIATE(std::complex<float>);
#  endif
#endif


#ifdef ENABLE_DOUBLE_PRECISION
INSTANTIATE(double);
#  ifdef ENABLE_COMPLEX_BASIS_FUNCTIONS
INSTANTIATE(std::complex<double>);
#  endif
#endif

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(Space);

} // namespace Bempp
