#ifndef bempp_boundary_operator_hpp
#define bempp_boundary_operator_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "transposition_mode.hpp"

#include <boost/mpl/set.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/utility/enable_if.hpp>
#include <string>

#include <tbb/mutex.h>

namespace Bempp
{

template <typename BasisFunctionType> class Space;
template <typename ResultType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class AbstractBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class Context;
template <typename BasisFunctionType, typename ResultType> class GridFunction;

template <typename BasisFunctionType, typename ResultType>
class BoundaryOperator
{
public:
    BoundaryOperator(const shared_ptr<const Context<
                     BasisFunctionType, ResultType> >& context,
                     const shared_ptr<const AbstractBoundaryOperator<
                     BasisFunctionType, ResultType> >& abstractOp);

    shared_ptr<const AbstractBoundaryOperator<BasisFunctionType, ResultType> >
    abstractOperator() const;

    shared_ptr<const Context<BasisFunctionType, ResultType> > context() const;

    shared_ptr<const DiscreteBoundaryOperator<ResultType> > weakForm() const;

    shared_ptr<const Space<BasisFunctionType> > domain() const;
    shared_ptr<const Space<BasisFunctionType> > range() const;
    shared_ptr<const Space<BasisFunctionType> > dualToRange() const;

    std::string label() const;

    /** \brief Set <tt>y_inout := alpha * A * x_in + beta * y_inout</tt>, where
     *  \c A is this operator. */
    void apply(const TranspositionMode trans,
               const GridFunction<BasisFunctionType, ResultType>& x_in,
               GridFunction<BasisFunctionType, ResultType>& y_inout,
               ResultType alpha, ResultType beta) const;

private:
    shared_ptr<const Context<BasisFunctionType, ResultType> > m_context;
    shared_ptr<const AbstractBoundaryOperator<BasisFunctionType, ResultType> >
    m_abstractOp;
    mutable shared_ptr<const DiscreteBoundaryOperator<ResultType> > m_weakForm;
    mutable tbb::mutex m_weakFormMutex;
};

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType> operator+(
        const BoundaryOperator<BasisFunctionType, ResultType>& op1,
        const BoundaryOperator<BasisFunctionType, ResultType>& op2);

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType> operator-(
        const BoundaryOperator<BasisFunctionType, ResultType>& op1,
        const BoundaryOperator<BasisFunctionType, ResultType>& op2);

// This type machinery is needed to disambiguate between this operator and
// the one taking a BoundaryOperator and a GridFunction
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    BoundaryOperator<BasisFunctionType, ResultType> >::type
operator*(
        const BoundaryOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar);

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
BoundaryOperator<BasisFunctionType, ResultType> operator*(
        const ScalarType& scalar,
        const BoundaryOperator<BasisFunctionType, ResultType>& op);

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
BoundaryOperator<BasisFunctionType, ResultType> operator/(
        const BoundaryOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar);

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const BoundaryOperator<BasisFunctionType, ResultType>& op,
        const GridFunction<BasisFunctionType, ResultType>& fun);

//template <typename BasisFunctionType, typename ResultType>
//BoundaryOperatorComposition<BasisFunctionType, ResultType> operator*(
//        const BoundaryOperator<BasisFunctionType, ResultType>& op1,
//        const BoundaryOperator<BasisFunctionType, ResultType>& op2);

} // namespace Bempp

#endif
