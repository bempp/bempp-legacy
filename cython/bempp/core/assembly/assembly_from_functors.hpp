#ifndef bempp_cython_assembly_assembly_from_functors_hpp
#define bempp_cython_assembly_assembly_from_functors_hpp

#include "bempp/assembly/general_elementary_local_operator.hpp"
#include "bempp/core/fiber/shape_transformation_functors.hpp"
#include "bempp/core/fiber/local_integrand_functors.hpp"

namespace Bempp {

    using namespace Fiber;

    inline shared_ptr<ElementaryLocalOperator<double, double>>
        abstract_local_operator_from_functors(
                const shared_ptr<const Space<double>> &domain,
                const shared_ptr<const Space<double>> &range,
                const shared_ptr<const Space<double>> &dualToRange,
                const std::string &label, int symmetry,
                const ShapeTransformationFunctorContainer& testFunctor,
                const ShapeTransformationFunctorContainer& trialFunctor,
                const LocalIntegrandFunctorContainer& integrandFunctor)
        {

        return shared_ptr<ElementaryLocalOperator<double, double>>(
                    new GeneralElementaryLocalOperator<double, double>(
                        domain, range, dualToRange, label, symmetry,
                        testFunctor, trialFunctor, integrandFunctor));

        }



}


#endif
