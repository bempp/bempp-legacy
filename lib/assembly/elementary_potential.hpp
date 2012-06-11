#ifndef bempp_elementary_potential_hpp
#define bempp_elementary_potential_hpp

#include "potential.hpp"

namespace Fiber
{

template <typename CoordinateType> class CollectionOfBasisTransformations;
template <typename KernelType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class KernelTrialIntegral;
template <typename ResultType> class EvaluatorForIntegralOperators;

} // namespace Bempp

namespace Bempp
{

template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
class ElementaryPotential : public Potential<BasisFunctionType_, ResultType_>
{
    typedef Potential<BasisFunctionType_, ResultType_> Base;
public:
    typedef typename Base::BasisFunctionType BasisFunctionType;
    typedef typename Base::ResultType ResultType;
    typedef KernelType_ KernelType;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::LocalAssemblerFactory LocalAssemblerFactory;
    typedef Fiber::EvaluatorForIntegralOperators<ResultType> Evaluator;
    typedef Fiber::CollectionOfBasisTransformations<CoordinateType>
    CollectionOfBasisTransformations;
    typedef Fiber::CollectionOfKernels<KernelType> CollectionOfKernels;
    typedef Fiber::KernelTrialIntegral<BasisFunctionType, KernelType, ResultType>
    KernelTrialIntegral;

    virtual std::auto_ptr<InterpolatedFunction<ResultType_> > evaluateOnGrid(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const Grid& evaluationGrid,
            const LocalAssemblerFactory& assemblerFactory,
            const EvaluationOptions& options) const;

    virtual arma::Mat<ResultType_> evaluateAtPoints(
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

    virtual const CollectionOfKernels& kernels() const = 0;
    virtual const CollectionOfBasisTransformations&
    trialTransformations() const = 0;
    virtual const KernelTrialIntegral& integral() const = 0;
};

} // namespace Bempp

#endif
