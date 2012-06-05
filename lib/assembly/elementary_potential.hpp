#ifndef bempp_elementary_potential_hpp
#define bempp_elementary_potential_hpp

#include "potential.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
class ElementaryPotential : public Potential<BasisFunctionType, ResultType>
{
public:
    virtual std::auto_ptr<InterpolatedFunction<ResultType> > applyOffSurface(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const Grid& evaluationGrid,
            const LocalAssemblerFactory& factory,
            const EvaluationOptions& options) const;

    virtual arma::Mat<ResultType> applyOffSurface(
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

#endif // bempp_elementary_potential_hpp
