
#ifndef bempp_cython_local_integrand_functors_hpp
#define bempp_cython_local_integrand_functors_hpp

#include "bempp/fiber/collection_of_3d_arrays.hpp"
#include "bempp/fiber/geometrical_data.hpp"
#include "bempp/fiber/simple_test_trial_integrand_functor.hpp"

namespace Fiber {



class LocalIntegrandFunctorBase
{
    public:

        virtual void addGeometricalDependencies(size_t &geomDeps) const = 0;

        inline virtual ~LocalIntegrandFunctorBase() {};

        virtual void evaluate(const ConstGeometricalDataSlice<double> &geomData,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& testValues,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& trialValues)
            const = 0;

};

template <typename Functor>
class ConcreteLocalIntegrandFunctor :
    public LocalIntegrandFunctorBase
{

    public:


        ConcreteLocalIntegrandFunctor(const Functor& functor) :
            m_functor(functor) {}

        virtual ~ConcreteLocalIntegrandFunctor() {};

        void addGeometricalDependencies(size_t &geomDeps) const override
        {
            m_functor.addGeometricalDependencies(geomDeps);
        }

        void evaluate(const ConstGeometricalDataSlice<double> &geomData,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& testValues,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& trialValues)
            const override {

                m_functor.evaluate(geomData, testValues, trialValues);
            } 

    private:

        Functor m_functor;

};

class LocalIntegrandFunctorContainer
{
    public:
        typedef double CoordinateType;
        typedef double ResultType;
        typedef double BasisFunctionType;

        inline LocalIntegrandFunctorContainer(const shared_ptr<const LocalIntegrandFunctorBase>&
                functor) : m_functor(functor) {}

        inline void addGeometricalDependencies(size_t &geomDeps) const {

            m_functor->addGeometricalDependencies(geomDeps);

        }

        inline void evaluate(const ConstGeometricalDataSlice<double> &geomData,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& testValues,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& trialValues)
            const {

                m_functor->evaluate(geomData, testValues, trialValues);
            } 


    private:

        shared_ptr<const LocalIntegrandFunctorBase> m_functor;

};



inline LocalIntegrandFunctorContainer* simpleTestTrialIntegrandFunctor(){

    return new LocalIntegrandFunctorContainer(
                shared_ptr<LocalIntegrandFunctorBase>(
                    new ConcreteLocalIntegrandFunctor<SimpleTestTrialIntegrandFunctor<double, double>>(
                        SimpleTestTrialIntegrandFunctor<double, double>())
                    )
                );

}



}


#endif
