#include "elementary_linear_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<typename ElementaryLinearOperator<BasisFunctionType, ResultType>::LocalAssembler>
ElementaryLinearOperator<BasisFunctionType, ResultType>::makeAssemblerFromScratch(
        const LocalAssemblerFactory& assemblerFactory,
        const AssemblyOptions& options) const
{
    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;

    shared_ptr<RawGridGeometry> testRawGeometry, trialRawGeometry;
    shared_ptr<GeometryFactory> testGeometryFactory, trialGeometryFactory;
    shared_ptr<Fiber::OpenClHandler> openClHandler;
    shared_ptr<BasisPtrVector> testBases, trialBases;
    bool cacheSingularIntegrals;

    collectDataForAssemblerConstruction(options,
                                        testRawGeometry, trialRawGeometry,
                                        testGeometryFactory, trialGeometryFactory,
                                        testBases, trialBases,
                                        openClHandler, cacheSingularIntegrals);

    return makeAssembler(assemblerFactory, testGeometryFactory, testRawGeometry,
                         // TODO: add parameters for trial*
                         testBases, trialBases, openClHandler,
                         options.parallelisationOptions(), cacheSingularIntegrals);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(ElementaryLinearOperator);

} // namespace Bempp
