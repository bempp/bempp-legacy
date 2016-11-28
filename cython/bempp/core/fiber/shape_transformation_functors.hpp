#ifndef bempp_cython_shape_transformation_functors_base_hpp
#define bempp_cython_shape_transformation_functors_base_hpp

#include "bempp/fiber/shared_ptr.hpp"
#include "bempp/fiber/basis_data.hpp"
#include "bempp/fiber/geometrical_data.hpp"
#include "bempp/fiber/collection_of_3d_arrays.hpp"
#include "bempp/fiber/shape_transformation_functor_wrappers.hpp"
#include "bempp/fiber/scalar_function_value_functor.hpp"
#include "bempp/fiber/surface_grad_3d_functor.hpp"
#include "bempp/fiber/surface_curl_3d_functor.hpp"
#include "bempp/fiber/surface_div_3d_functor.hpp"
#include "bempp/fiber/hcurl_function_value_functor.hpp"
#include "bempp/fiber/hdiv_function_value_functor.hpp"
#include "bempp/fiber/hcurl_surface_curl_functor.hpp"
#include "bempp/fiber/scalar_function_value_times_normal_functor.hpp"

namespace Fiber {


class ShapeTransformationFunctorBase
{
    public:

        inline virtual ~ShapeTransformationFunctorBase() {};

        virtual size_t transformationCount() const = 0;
        
        virtual int argumentDimension() const = 0; 

        virtual int resultDimension(size_t transformationIndex) const = 0;

        virtual void addDependencies(size_t &basisDeps, size_t &geomDeps) const = 0; 

        virtual void evaluate(const ConstBasisDataSlice<double> &basisData,
                      const ConstGeometricalDataSlice<double> &geomData,
                      CollectionOf1dSlicesOf3dArrays<double> &result) const = 0;

};

template <typename Functor>
class ConcreteShapeTransformationFunctor :
    public ShapeTransformationFunctorBase
{

    public:

        typedef double CoordinateType;

        ConcreteShapeTransformationFunctor(const Functor& functor) :
            m_functor(functor) {}

        inline virtual ~ConcreteShapeTransformationFunctor() {};

        size_t transformationCount() const override { 
            return m_functor.transformationCount();
        }

        int argumentDimension() const override {

            return m_functor.argumentDimension();
        }

        int resultDimension(size_t transformationIndex) const override {

            return m_functor.resultDimension(transformationIndex);

        }

        void addDependencies(size_t& basisDeps, size_t &geomDeps) const override {
            m_functor.addDependencies(basisDeps, geomDeps);
        }

        void evaluate(const ConstBasisDataSlice<double> &basisData,
                const ConstGeometricalDataSlice<double> &geomData,
                CollectionOf1dSlicesOf3dArrays<double> & result) const override
        {
            m_functor.evaluate(basisData, geomData, result);
        }

    private:

        Functor m_functor;

};

class ShapeTransformationFunctorContainer
{
    public:
        typedef double CoordinateType;

        inline ShapeTransformationFunctorContainer(const shared_ptr<const ShapeTransformationFunctorBase>&
                functor) : m_functor(functor) {}

        inline size_t transformationCount() const { return m_functor->transformationCount();}
        
        inline int argumentDimension() const {return m_functor->argumentDimension();}

        inline int resultDimension(size_t transformationIndex) const {
            return m_functor->resultDimension(transformationIndex);
        }

        inline void addDependencies(size_t &basisDeps, size_t &geomDeps) const {
            m_functor->addDependencies(basisDeps, geomDeps);

        }

        inline void evaluate(const ConstBasisDataSlice<double> &basisData,
                      const ConstGeometricalDataSlice<double> &geomData,
                      CollectionOf1dSlicesOf3dArrays<double> &result) const {

            m_functor->evaluate(basisData, geomData, result);

        }

        inline void evaluate(const ConstBasisDataSlice<std::complex<double>> &basisData,
                      const ConstGeometricalDataSlice<double> &geomData,
                      CollectionOf1dSlicesOf3dArrays<std::complex<double>> &result) const {
            // Only necessary because CollectionOfShapeTransformations always wants to
            // implement evaluate for complex types.

            throw std::runtime_error("ShapeTransformationFunctorContainer:evaluate(): "
                                     "Not implemented for complex types.");

        }

       

    private:

        shared_ptr<const ShapeTransformationFunctorBase> m_functor;

};

inline ShapeTransformationFunctorContainer* scalarFunctionValueFunctor(){

    return new ShapeTransformationFunctorContainer(
                shared_ptr<ShapeTransformationFunctorBase>(
                    new ConcreteShapeTransformationFunctor<ScalarFunctionValueFunctor<double>>(
                        ScalarFunctionValueFunctor<double>())
                    )
                );

}

inline ShapeTransformationFunctorContainer* surfaceGrad3dFunctor(){

    return new ShapeTransformationFunctorContainer(
                shared_ptr<ShapeTransformationFunctorBase>(
                    new ConcreteShapeTransformationFunctor<SurfaceGrad3dFunctor<double>>(
                        SurfaceGrad3dFunctor<double>())
                    )
                );

}

inline ShapeTransformationFunctorContainer* surfaceDiv3dFunctor(){

    return new ShapeTransformationFunctorContainer(
                shared_ptr<ShapeTransformationFunctorBase>(
                    new ConcreteShapeTransformationFunctor<SurfaceDiv3dFunctor<double>>(
                        SurfaceDiv3dFunctor<double>())
                    )
                );

}

inline ShapeTransformationFunctorContainer* surfaceCurl3dFunctor(){

    return new ShapeTransformationFunctorContainer(
                shared_ptr<ShapeTransformationFunctorBase>(
                    new ConcreteShapeTransformationFunctor<SurfaceCurl3dFunctor<double>>(
                        SurfaceCurl3dFunctor<double>())
                    )
                );

}


inline ShapeTransformationFunctorContainer* hcurlFunctionValueFunctor(){

    return new ShapeTransformationFunctorContainer(
                shared_ptr<ShapeTransformationFunctorBase>(
                    new ConcreteShapeTransformationFunctor<HcurlFunctionValueFunctor<double>>(
                        HcurlFunctionValueFunctor<double>())
                    )
                );

}

inline ShapeTransformationFunctorContainer* hdivFunctionValueFunctor(){

    return new ShapeTransformationFunctorContainer(
                shared_ptr<ShapeTransformationFunctorBase>(
                    new ConcreteShapeTransformationFunctor<HdivFunctionValueFunctor<double>>(
                        HdivFunctionValueFunctor<double>())
                    )
                );

}

inline ShapeTransformationFunctorContainer* hcurlSurfaceCurlFunctor(){

    return new ShapeTransformationFunctorContainer(
                shared_ptr<ShapeTransformationFunctorBase>(
                    new ConcreteShapeTransformationFunctor<HcurlSurfaceCurlFunctor<double>>(
                        HcurlSurfaceCurlFunctor<double>())
                    )
                );

}

inline ShapeTransformationFunctorContainer* scalarFunctionValueTimesNormalFunctor(){

    return new ShapeTransformationFunctorContainer(
                shared_ptr<ShapeTransformationFunctorBase>(
                    new ConcreteShapeTransformationFunctor<ScalarFunctionValueTimesNormalFunctor<double>>(
                        ScalarFunctionValueTimesNormalFunctor<double>())
                    )
                );

}

}

#endif
