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

#include "bempp/common/config_trilinos.hpp"
#include "discrete_fmm_boundary_operator.hpp"
#include "octree.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/type_traits/is_complex.hpp>


#ifdef WITH_TRILINOS
#include <Thyra_DefaultSpmdVectorSpace_decl.hpp>
#endif

namespace Bempp
{


template <typename ValueType>
DiscreteFmmBoundaryOperator<ValueType>::
DiscreteFmmBoundaryOperator(
		unsigned int rowCount, unsigned int columnCount,
		const shared_ptr<Octree<ValueType> > &octree,
		int symmetry_) :
#ifdef WITH_TRILINOS
    m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(columnCount)),
    m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(rowCount)),
#else
    m_rowCount(rowCount), m_columnCount(columnCount),
#endif
	m_octree(octree),
	m_symmetry(symmetry_)
{
}

template <typename ValueType>
arma::Mat<ValueType>
DiscreteFmmBoundaryOperator<ValueType>::
asMatrix() const
{
}

template <typename ValueType>
unsigned int
DiscreteFmmBoundaryOperator<ValueType>::
rowCount() const
{
#ifdef WITH_TRILINOS
    return m_rangeSpace->dim();
#else
    return m_rowCount;
#endif
}

template <typename ValueType>
unsigned int
DiscreteFmmBoundaryOperator<ValueType>::
columnCount() const
{
#ifdef WITH_TRILINOS
    return m_domainSpace->dim();
#else
    return m_columnCount;
#endif
}

template <typename ValueType>
void
DiscreteFmmBoundaryOperator<ValueType>::
addBlock(const std::vector<int>& rows,
         const std::vector<int>& cols,
         const ValueType alpha,
         arma::Mat<ValueType>& block) const
{
    throw std::runtime_error("DiscreteFmmBoundaryOperator::"
                             "addBlock(): not implemented yet");
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> >
DiscreteFmmBoundaryOperator<ValueType>::asDiscreteFmmBoundaryOperator(
        double eps, int maximumRank, bool interleave) const
{
    return this->shared_from_this(); // this-> needed for template name resolution.
}

template <typename ValueType>
const DiscreteFmmBoundaryOperator<ValueType>&
DiscreteFmmBoundaryOperator<ValueType>::castToFmm(
        const DiscreteBoundaryOperator<ValueType>& discreteOperator)
{
    return dynamic_cast<const DiscreteFmmBoundaryOperator<ValueType>&>(
                discreteOperator);
}

template <typename ValueType>
shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> >
DiscreteFmmBoundaryOperator<ValueType>::castToFmm(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >&
        discreteOperator)
{
    shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> > result =
        boost::dynamic_pointer_cast<const DiscreteFmmBoundaryOperator<ValueType> >(
                discreteOperator);
    if (result.get()==0 && (discreteOperator.get()!=0)) throw std::bad_cast();
    return result;
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteFmmBoundaryOperator<ValueType>::
domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteFmmBoundaryOperator<ValueType>::
range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool
DiscreteFmmBoundaryOperator<ValueType>::
opSupportedImpl(Thyra::EOpTransp M_trans) const
{
    // TODO: implement remaining variant (conjugate)
    return (M_trans == Thyra::NOTRANS || M_trans == Thyra::TRANS ||
            M_trans == Thyra::CONJTRANS);
}
#endif // WITH_TRILINOS

template <typename ValueType>
void
DiscreteFmmBoundaryOperator<ValueType>::
applyBuiltInImpl(const TranspositionMode trans,
                 const arma::Col<ValueType>& x_in,
                 arma::Col<ValueType>& y_inout,
                 const ValueType alpha,
                 const ValueType beta) const
{
/*
    if (!transposed)
        m_rangePermutation.unpermuteVector(permutedResult, y_inout);
    else
        m_domainPermutation.unpermuteVector(permutedResult, y_inout);
*/
    if (trans != NO_TRANSPOSE && trans != TRANSPOSE && trans != CONJUGATE_TRANSPOSE)
        throw std::runtime_error(
                "DiscreteFmmBoundaryOperator::applyBuiltInImpl(): "
                "transposition modes other than NO_TRANSPOSE and "
                "CONJUGATE_TRANSPOSE are not supported");

	bool transposed = (trans & TRANSPOSE);

    if ((!transposed && (columnCount() != x_in.n_rows ||
                         rowCount() != y_inout.n_rows)) ||
            (transposed && (rowCount() != x_in.n_rows ||
                            columnCount() != y_inout.n_rows)))
        throw std::invalid_argument(
                "DiscreteFmmBoundaryOperator::applyBuiltInImpl(): "
                "incorrect vector length");

	arma::Col<ValueType> y_in = y_inout;
	y_inout.fill(0.0);
	m_octree->matrixVectorProduct(x_in, y_inout, trans);
	y_inout = alpha*y_inout + beta*y_in;
}

template <typename ValueType>
double
DiscreteFmmBoundaryOperator<ValueType>::eps() const
{
    return m_eps;
}

template <typename ValueType>
int
DiscreteFmmBoundaryOperator<ValueType>::maximumRank() const
{
    return m_maximumRank;
}

template <typename ValueType>
int
DiscreteFmmBoundaryOperator<ValueType>::symmetry() const
{
    return m_symmetry;
}

template <typename ValueType>
const IndexPermutation&
DiscreteFmmBoundaryOperator<ValueType>::domainPermutation() const
{
    //return m_domainPermutation;
}

template <typename ValueType>
const IndexPermutation&
DiscreteFmmBoundaryOperator<ValueType>::rangePermutation() const
{
    //return m_rangePermutation;
}

template <typename ValueType>
const ParallelizationOptions&
DiscreteFmmBoundaryOperator<ValueType>::parallelizationOptions() const
{
    return m_parallelizationOptions;
}

// Global routines

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > fmmOperatorSum(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
        double eps, int maximumRank)
{

}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator(
        const ValueType& multiplier,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        const ValueType& multiplier)
{
    return scaledFmmOperator(multiplier, op);
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteFmmBoundaryOperator);

#define INSTANTIATE_FREE_FUNCTIONS(RESULT) \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT> > \
        fmmOperatorSum( \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op1, \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op2, \
            double eps, int maximumRank); \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT> > \
        scaledFmmOperator( \
            const RESULT& multiplier, \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op); \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT> > \
        scaledFmmOperator( \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op, \
            const RESULT& multiplier)

#if defined(ENABLE_SINGLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(float);
#endif

#if defined(ENABLE_SINGLE_PRECISION) && (defined(ENABLE_COMPLEX_BASIS_FUNCTIONS) || defined(ENABLE_COMPLEX_KERNELS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<float>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(double);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && (defined(ENABLE_COMPLEX_BASIS_FUNCTIONS) || defined(ENABLE_COMPLEX_KERNELS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<double>);
#endif

} // namespace Bempp
