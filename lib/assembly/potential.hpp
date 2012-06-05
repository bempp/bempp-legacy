#ifndef bempp_potential_hpp
#define bempp_potential_hpp

#include "../fiber/local_assembler_factory.hpp"

#include <armadillo>
#include <memory>

namespace Bempp
{

class EvaluationOptions;
class Grid;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename ResultType> class InterpolatedFunction;

template <typename BasisFunctionType, typename ResultType>
class Potential
{
public:
    typedef Fiber::LocalAssemblerFactory<BasisFunctionType, ResultType, GeometryFactory>
    LocalAssemblerFactory;

    virtual std::auto_ptr<InterpolatedFunction<ResultType> > applyOffSurface(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const Grid& evaluationGrid,
            const LocalAssemblerFactory& factory,
            const EvaluationOptions& options) const = 0;

    virtual arma::Mat<ResultType> applyOffSurface(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const arma::Mat<CoordinateType>& evaluationPoints,
            const LocalAssemblerFactory& assemblerFactory,
            const EvaluationOptions& options) const = 0;
};

} // namespace Bempp

#endif
