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


#ifndef bempp_identity_operator_hpp
#define bempp_identity_operator_hpp

#include "../common/common.hpp"

#include "elementary_abstract_boundary_operator.hpp"

#include "abstract_boundary_operator_id.hpp"
#include <boost/scoped_ptr.hpp>

namespace Fiber
{

template <typename ResultType> class LocalAssemblerForOperators;

} // namespace Fiber

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType> class BoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class IdentityOperator;

template <typename BasisFunctionType, typename ResultType>
class IdentityOperatorId : public AbstractBoundaryOperatorId
{
public:
    IdentityOperatorId(const IdentityOperator<BasisFunctionType, ResultType>& op);
    virtual size_t hash() const;
    virtual bool isEqual(const AbstractBoundaryOperatorId &other) const;

private:
    const Space<BasisFunctionType>* m_domain;
    const Space<BasisFunctionType>* m_range;
    const Space<BasisFunctionType>* m_dualToRange;
};

/** \ingroup identity
 *  \brief Identity operator.
 *
 *  Let \f$X\f$ and \f$Y\f$ be two function spaces defined on the same grid. If
 *  \f$X \supset Y\f$, an instance of IdentityOperator with domain \f$X\f$
 *  and range \f$Y\f$ represents the orthogonal projection operator from
 *  \f$X\f$ to \f$Y\f$. If \f$X \subset Y\f$, it represents the inclusion
 *  operator from \f$X\f$ to \f$Y\f$. In the special case of \f$X = Y\f$, we
 *  have the standard identity operator. In BEM++ we (ab)use the term "identity
 *  operator" to refer to all these three cases.
 *
 *  See AbstractBoundaryOperator for the documentation of the template
 *  parameters.
 *
 *  Use identityOperator() to create a BoundaryOperator object wrapping
 *  an identity operator.
 */
template <typename BasisFunctionType_, typename ResultType_>
class IdentityOperator :
        public ElementaryAbstractBoundaryOperator<BasisFunctionType_, ResultType_>
{
    typedef ElementaryAbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc ElementaryAbstractBoundaryOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc ElementaryAbstractBoundaryOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc ElementaryAbstractBoundaryOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc ElementaryAbstractBoundaryOperator::QuadratureStrategy */
    typedef typename Base::QuadratureStrategy QuadratureStrategy;
    /** \copydoc ElementaryAbstractBoundaryOperator::LocalAssembler */
    typedef typename Base::LocalAssembler LocalAssembler;

    /** \brief Constructor.
     *
     *  \param[in] domain
     *    Function space being the domain of the operator.
     *  \param[in] range
     *    Function space being the range of the operator.
     *  \param[in] dualToRange
     *    Function space dual to the the range of the operator.
     *  \param[in] label
     *    Textual label of the operator. If empty, a unique label is generated
     *    automatically.
     *  \param[in] symmetry
     *    Symmetry of the weak form of the operator. Can be any combination of
     *    the flags defined in the enumeration type Symmetry.
     *    If set to AUTO_SYMMETRY (default), the symmetry is determined
     *    automatically by checking whether its domain and space dual to its
     *    range are equal. If so, the operator is marked as Hermitian,
     *    and if the basis functions are real-valued, also as symmetric.
     *
     *  All the three spaces must be defined on the same grid. */
    IdentityOperator(const shared_ptr<const Space<BasisFunctionType> >& domain,
                     const shared_ptr<const Space<BasisFunctionType> >& range,
                     const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                     const std::string& label = "",
                     int symmetry = AUTO_SYMMETRY);
    IdentityOperator(const IdentityOperator& other);
    virtual ~IdentityOperator();

    /** \brief Return the identifier of this operator.
     *
     *  Identity operators are treated as equivalent if they have the same domain,
     *  range and dual to range. */
    virtual shared_ptr<const AbstractBoundaryOperatorId> id() const;

    /** \brief Return true. */
    virtual bool isLocal() const;

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

private:
    virtual std::auto_ptr<LocalAssembler> makeAssemblerImpl(
            const QuadratureStrategy& quadStrategy,
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
            const ParallelizationOptions& parallelizationOptions,
            bool cacheSingularIntegrals) const;

    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInternalImpl(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInDenseMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInSparseMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

private:
    struct Impl;
    boost::scoped_ptr<Impl> m_impl;
    shared_ptr<const AbstractBoundaryOperatorId> m_id;
};

/** \relates IdentityOperator
 *  \brief Construct a BoundaryOperator object wrapping an IdentityOperator.
 *
 *  This convenience function constructs an abstract identity operator and wraps
 *  it in a BoundaryOperator object.
 *
 *  \param[in] context
 *    A Context object that will be used to build the weak form of the
 *    identity operator when necessary.
 *  \param[in] domain
 *    Function space being the domain of the identity operator.
 *  \param[in] range
 *    Function space being the range of the identity operator.
 *  \param[in] dualToRange
 *    Function space dual to the the range of the identity operator.
 *  \param[in] label
 *    Textual label of the identity operator (optional, used for debugging).
 *  \param[in] symmetry
 *    Symmetry of the weak form of the operator. Can be any combination of
 *    the flags defined in the enumeration type Symmetry.
 *    If set to AUTO_SYMMETRY (default), the symmetry is determined
 *    automatically by checking whether its domain and space dual to its
 *    range are equal. If so, the operator is marked as Hermitian.
 */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
identityOperator(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& domain,
                 const shared_ptr<const Space<BasisFunctionType> >& range,
                 const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                 const std::string& label = "",
                 int symmetry = AUTO_SYMMETRY);

} // namespace Bempp

#endif
