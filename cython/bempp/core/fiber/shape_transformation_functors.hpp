#ifndef bempp_cython_shape_transformation_functors_base_hpp
#define bempp_cython_shape_transformation_functors_base_hpp

#include "bempp/fiber/shared_ptr.hpp"
#include "bempp/fiber/basis_data.hpp"
#include "bempp/fiber/geometrical_data.hpp"
#include "bempp/fiber/collection_of_3d_arrays.hpp"
#include "bempp/fiber/shape_transformation_functor_wrappers.hpp"
#include "bempp/fiber/scalar_function_value_functor.hpp"

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

        void evaluate(const ConstBasisDataSlice<double> &basisData,
                      const ConstGeometricalDataSlice<double> &geomData,
                      CollectionOf1dSlicesOf3dArrays<double> &result) const {

            m_functor->evaluate(basisData, geomData, result);

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

}

#endif
