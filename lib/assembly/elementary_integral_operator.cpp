#include "elementary_integral_operator.hpp"

#include "assembly_options.hpp"
#include "kernel_adapter.hpp"
#include "quadrature_selector.hpp"

#include "../common/multidimensional_arrays.hpp"
#include "../common/not_implemented_error.hpp"
#include "../fiber/double_integrator.hpp"
#include "../grid/entity_pointer.hpp"
#include "../grid/entity.hpp"
#include "../grid/geometry_adapter.hpp"
#include "../space/function_family_adapter.hpp"
#include "../space/space.hpp"

#include <armadillo>
// #include <boost/unordered_set.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <set>

namespace Bempp
{

template <typename ValueType>
void ElementaryIntegralOperator<ValueType>::evaluateLocalWeakForms(
                const std::vector<const EntityPointer<0>*>& testElements,
                const std::vector<const EntityPointer<0>*>& trialElements,
                const Space<ValueType>& testSpace,
                const Space<ValueType>& trialSpace,
                const QuadratureSelector<ValueType>& quadSelector,
                const AssemblyOptions& options,
                Array2D<arma::Mat<ValueType> >& result) const
{
    // TODO: if the operator, spaces or geometry do not support OpenCL,
    // detect it (via catching an exception), warn user and fall back on
    // the non-OpenCL implementation
    if (options.useOpenCl)
        evaluateLocalWeakFormsWithOpenCl(testElements, trialElements,
                                         testSpace, trialSpace,
                                         quadSelector, options, result);
    else
        evaluateLocalWeakFormsWithoutOpenCl(testElements, trialElements,
                                            testSpace, trialSpace,
                                            quadSelector, options, result);
}

template <typename ValueType>
void ElementaryIntegralOperator<ValueType>::evaluateLocalWeakFormsWithOpenCl(
        const std::vector<const EntityPointer<0>*>& testElements,
        const std::vector<const EntityPointer<0>*>& trialElements,
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const QuadratureSelector<ValueType>& quadSelector,
        const AssemblyOptions& options,
        Array2D<arma::Mat<ValueType> >& result) const
{
    throw NotImplementedError("ElementaryIntegralOperator::"
                              "evaluateLocalWeakFormsWithOpenCl(): "
                              "not implemented yet");

    // TODO: update this comment in view of the Fiber interface

    // Needed functions:
    //
    // - Space:
    //   void evaluateBasisFunction(int elementType, const Point* local,
    //     int localDofIndex, Float* result);
    //   void evaluateBasisFunctionDerivative(int elementType, const Point* local,
    //     int localDofIndex, int direction, Float* result);
    //   int basisFunctionCount(int elementType); // maybe unnecessary
    //   void evaluateShapeFunction(int elementType, const Point* local,
    //     Float basisFunction, Float basisFunctionDerivativeDir0, Float basisFunctionDerivativeDir1,
    //     Float* jacobianMatrix, Float* result); // jacobianMatrix might be null if not needed
    //   bool needsJacobianMatrix(); // for shape function evaluation
    //   bool needsBasisFunctionDerivatives(); // for shape function evaluation
    //
    //
    // - Geometry:
    //   void getGeometricInformation(int elementType, Float* geometryData,
    //     const Point* local,
    //     Float* integrationElement, Float* global, Float* jacobianMatrix);
    //     // jacobianMatrix can be null -> is then not evaluated
    //   int geometryDataLength(); // number of elements of geometryData that
    //     // should be supplied to getGeometricInformation
    //
    // - IntegralOperator:
    //   bool needsJacobianMatrix();
    //   void evaluateKernel(const Point* testPoint, const Point* quadPoint,
    //     const Float* jacobianMatrix, Float* result); // result might be a tensor!
    //   void evaluateOperatorIntegrand(const Point* testPoint, const Point* quadPoint,
    //     Float basisFunction, Float basisFunctionDerivativeDir0,
    //     Float basisFunctionDerivativeDir1, const Float* jacobianMatrix,
    //     Float integrationElement, Float* result);
    //   void evaluateWeakFormIntegrand(const PointPair* quadPoint,
    //     Float basisFunctionAtTestPoint, Float basisFunctionDerivativeDir0,
    //     Float basisFunctionDerivativeDir1, const Float* jacobianMatrix,
    //     Float integrationElement, Float* result);
    //
    // - the Geometry class should provide a function
    //     void getOpaqueInfo(arma::Col<ValueType>& info);
    //   that will fill info with the data to be passed to OpenCL and from which
    //   (and possibly the geometryType) the geometry can be reconstructed.
    //   E.g. for a flat triangular element those might be the 3D coordinates
    //   of the three vertices.
    //   For a curved element -- some additional information will be needed.
}

template <typename ValueType>
void ElementaryIntegralOperator<ValueType>::evaluateLocalWeakFormsWithoutOpenCl(
        const std::vector<const EntityPointer<0>*>& testElements,
        const std::vector<const EntityPointer<0>*>& trialElements,
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const QuadratureSelector<ValueType>& quadSelector,
        const AssemblyOptions& options,
        Array2D<arma::Mat<ValueType> >& result) const
{
    using std::vector;

    const int testElementCount = testElements.size();
    const int trialElementCount = trialElements.size();

    Array2D<Fiber::QuadratureRule> quadRules;

    KernelAdapter<ValueType> kernelAdapter(kernel());

    // Construct adapters for all test geometries
    boost::ptr_vector<GeometryAdapter> owner;
    owner.reserve(testElementCount);
    for (int i = 0; i < testElementCount; ++i)
        owner.push_back(new GeometryAdapter(testElements[i]->entity().geometry()));
    vector<const GeometryAdapter*> testGeometries;
    testGeometries.reserve(testElementCount);
    for (int i = 0; i < testElementCount; ++i)
        testGeometries.push_back(&owner[i]);

    vector<const GeometryAdapter*> activeTrialGeometries(1);
    vector<const GeometryAdapter*> activeTestGeometries;
    activeTestGeometries.reserve(testElementCount);

    result.set_size(testElementCount, trialElementCount);

    vector<const EntityPointer<0>*> activeTrialElements(1);
    const EntityPointer<0>*& activeTrialElement = activeTrialElements[0];

    // Loop over trial elements
    for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
    {
        // Retrieve trial element variant and construct an adapter for its
        // geometry
        activeTrialElement = trialElements[trialIndex];
        ElementVariant trialElementVariant =
                trialSpace.elementVariant(activeTrialElement->entity());
        GeometryAdapter trialGeometry(activeTrialElement->entity().geometry());
        activeTrialGeometries[0] = &trialGeometry;

        // Choose appropriate quadrature rules for the integrals over the
        // current trial element and all test elements
        quadSelector.selectDoubleQuadratureRules(
                    testGeometries, activeTrialGeometries, quadRules);

        // Integration will proceed in batches of test elements having the same
        // "quadrature variant", i.e. quadrature rule and element variant (the
        // latter is important because it guarantees that all the elements have
        // the same basis functions)

        // First, find all the unique quadrature variants present
        typedef std::pair<Fiber::QuadratureRule, ElementVariant> QuadVariant;
        vector<QuadVariant> quadVariants; // Temporary vector
        quadVariants.reserve(testElementCount);
        for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
        {
            ElementVariant testElementVariant =
                    trialSpace.elementVariant(testElements[testIndex]->entity());
            quadVariants[testIndex] = QuadVariant(quadRules(testIndex, 0),
                                                  testElementVariant);
        }
        typedef std::set<QuadVariant> QuadVariantSet;
        // Set of unique quadrature variants
        QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());
        quadVariants.clear(); // We don't need the vector any more

        // Now loop over unique quadrature variants
        for (QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
             it != uniqueQuadVariants.end(); ++it)
        {
            Fiber::QuadratureRule quadRule = it->first;
            ElementVariant testElementVariant = it->second;

            // Get quadrature points and weights
            arma::Mat<ctype> testQuadPoints;
            arma::Mat<ctype> trialQuadPoints;
            vector<ValueType> quadWeights;
            quadSelector.doubleQuadratureRulePointsAndWeights(
                        quadRule, testQuadPoints, trialQuadPoints, quadWeights);
            int pointCount = quadWeights.size();

            // Find all the test elements for which quadrature should proceed
            // according to the current quadrature variant
            activeTestGeometries.clear();
            for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
                if (quadRules(testIndex, 0) == quadRule)
                    activeTestGeometries.push_back(testGeometries[testIndex]);

            // Set up test and trial function families
            std::auto_ptr<FunctionFamily<ValueType> > testFamily =
                    testFunctionFamily(testSpace, testElementVariant);
            std::auto_ptr<FunctionFamily<ValueType> > trialFamily =
                    trialFunctionFamily(trialSpace, trialElementVariant);
            FunctionFamilyAdapter<ValueType> testFamilyAdapter(*testFamily);
            FunctionFamilyAdapter<ValueType> trialFamilyAdapter(*trialFamily);

            // Integrate!
            // (Note: we pass a 3D rather than a 4D array as the last argument
            // to integrate() because we pass only one trial element.)
            arma::Cube<ValueType> localResult;
            Fiber::DoubleIntegrator::integrate(
                        pointCount,
                        testQuadPoints.begin(),
                        trialQuadPoints.begin(),
                        &*quadWeights.begin(),
                        activeTestGeometries.size(),
                        &*activeTestGeometries.begin(),
                        activeTrialGeometries.size(),
                        &*activeTrialGeometries.begin(),
                        testFamilyAdapter,
                        trialFamilyAdapter,
                        kernelAdapter,
                        localResult.begin());

            // Distribute the just calculated integrals into the result array
            // that will be returned to caller
            int i = 0;
            for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
                if (quadRules(testIndex,0) == quadRule)
                    result(testIndex, trialIndex) =
                            localResult(arma::span::all, arma::span(i++), arma::span::all);
        }
    }
}

#ifdef COMPILE_FOR_FLOAT
template class ElementaryIntegralOperator<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class ElementaryIntegralOperator<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class ElementaryIntegralOperator<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class ElementaryIntegralOperator<std::complex<double> >;
#endif

} // namespace Bempp
