#ifndef bempp_elementary_potential_hpp
#define bempp_elementary_potential_hpp

#include "potential.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename KernelType, typename ResultType>
class ElementaryPotential : public Potential<BasisFunctionType, ResultType>
{
    typedef Potential<BasisFunctionType, ResultType> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::LocalAssemblerFactory LocalAssemblerFactory;
    typedef Fiber::EvaluatorForIntegralOperators<ResultType> Evaluator;

    virtual std::auto_ptr<InterpolatedFunction<ResultType> > evaluateOnGrid(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const Grid& evaluationGrid,
            const LocalAssemblerFactory& assemblerFactory,
            const EvaluationOptions& options) const;

    virtual arma::Mat<ResultType> evaluateAtPoints(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const arma::Mat<CoordinateType>& evaluationPoints,
            const LocalAssemblerFactory& assemblerFactory,
            const EvaluationOptions& options) const;

private:
    std::auto_ptr<Evaluator>
    makeEvaluator(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const LocalAssemblerFactory& factory,
            const EvaluationOptions& options) const;

    virtual const Fiber::Kernel<KernelType>& kernel() const = 0;
    virtual const Fiber::ExpressionList<ResultType>& trialExpressionList() const = 0;
};

} // namespace Bempp

#endif
