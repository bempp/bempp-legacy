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


#ifndef bempp_null_operator_hpp
#define bempp_null_operator_hpp

#include "../common/common.hpp"

#include "abstract_boundary_operator.hpp"

#include "abstract_boundary_operator_id.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename BasisFunctionType, typename ResultType> class BoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class NullOperator;
/** \endcond */

template <typename BasisFunctionType, typename ResultType>
class BEMPP_DEPRECATED NullOperatorId : public AbstractBoundaryOperatorId
{
public:
    NullOperatorId(const NullOperator<BasisFunctionType, ResultType>& op);
    virtual size_t hash() const;
    virtual bool isEqual(const AbstractBoundaryOperatorId &other) const;

private:
    const Space<BasisFunctionType>* m_domain;
    const Space<BasisFunctionType>* m_range;
    const Space<BasisFunctionType>* m_dualToRange;
};

/** \ingroup null
 *  \brief Null operator.
 *
 *  This class represents an operator that always produces a zero function.
 *
 *  See AbstractBoundaryOperator for the documentation of the template
 *  parameters.
 *
 *  Use nullOperator() to create a BoundaryOperator object wrapping
 *  a null operator.
 */
template <typename BasisFunctionType_, typename ResultType_>
class NullOperator :
        public AbstractBoundaryOperator<BasisFunctionType_, ResultType_>
{
    typedef AbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc AbstractBoundaryOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc AbstractBoundaryOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc AbstractBoundaryOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc AbstractBoundaryOperator::QuadratureStrategy */
    typedef typename Base::QuadratureStrategy QuadratureStrategy;

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
     *    range are equal. If so, the operator is marked as symmetric and
     *    Hermitian. */
    NullOperator(const shared_ptr<const Space<BasisFunctionType> >& domain,
                     const shared_ptr<const Space<BasisFunctionType> >& range,
                     const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                     const std::string& label = "",
                     int symmetry = AUTO_SYMMETRY);
    NullOperator(const NullOperator& other);
    virtual ~NullOperator();

    /** \brief Return the identifier of this operator.
     *
     *  Null operators are treated as equivalent if they have the same domain,
     *  range and dual to range.
     *
     *  \deprecated This function is deprecated and will be removed in a future
     *  version of BEM++. */
    BEMPP_DEPRECATED virtual shared_ptr<const AbstractBoundaryOperatorId> id() const;

    /** \brief Return true. */
    virtual bool isLocal() const;

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

    shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    reallyAssembleWeakForm() const;

private:
    shared_ptr<const AbstractBoundaryOperatorId> m_id;
};

/** \relates NullOperator
 *  \brief Construct a BoundaryOperator object wrapping a NullOperator.
 *
 *  This convenience function constructs an abstract null operator and wraps
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
 *    Textual label of the operator (optional, used for debugging).
 *  \param[in] symmetry
 *    Symmetry of the weak form of the operator. Can be any combination of
 *    the flags defined in the enumeration type Symmetry.
 *    If set to AUTO_SYMMETRY (default), the symmetry is determined
 *    automatically by checking whether its domain and space dual to its
 *    range are equal. If so, the operator is marked as symetric and Hermitian.
 */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
nullOperator(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
             const shared_ptr<const Space<BasisFunctionType> >& domain,
             const shared_ptr<const Space<BasisFunctionType> >& range,
             const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
             const std::string& label = "",
             int symmetry = AUTO_SYMMETRY);

} // namespace Bempp

#endif
