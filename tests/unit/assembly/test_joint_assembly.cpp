// Copyright (C) 2011 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "bempp/common/config_ahmed.hpp"

#include "../type_template.hpp"
#include "../check_arrays_are_close.hpp"

#include "assembly/context.hpp"
#include "assembly/discrete_boundary_operator.hpp"
#include "assembly/identity_operator.hpp"
#include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "grid/grid_factory.hpp"
#include "space/piecewise_constant_scalar_space.hpp"
#include "space/piecewise_linear_continuous_scalar_space.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/type_traits/is_complex.hpp>
#include "grid/grid.hpp"

using namespace Bempp;

// Tests

namespace
{

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_V_plus_K
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
            context, domain, range, dualToRange) +
                laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_V_plus_3_K
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                3. * laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_3_V_plus_K
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return 3. * laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_3_V_plus_4_K
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return 3. * laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                4. * laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_V_plus_K_plus_Kp
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                laplace3dAdjointDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_3_V_plus_K_plus_Kp
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return 3. * laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                laplace3dAdjointDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_V_plus_3_K_plus_Kp
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                3. * laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                laplace3dAdjointDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_V
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_3_V
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return 3. * laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_I
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return identityOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_3_I
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return 3. * identityOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_I_plus_V
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return identityOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_V_plus_I
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                identityOperator<BasisFunctionType, ResultType>(
                                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_3_I_plus_V
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return 3. * identityOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_I_plus_3_V
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return identityOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                3. * laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_2_I_plus_V_plus_I_plus_4_K
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return 2. * identityOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                identityOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange) +
                2. * laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange);
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_K_plus_adj_Kp
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange, "DLP") +
                adjoint(laplace3dAdjointDoubleLayerBoundaryOperator<
                        BasisFunctionType, ResultType>(
                            context, dualToRange, dualToRange /* hack */, domain, "ADLP"));
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_K_plus_I_plus_4_adj_Kp
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange, "DLP") +
                identityOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange, "Id") +
                4. * adjoint(laplace3dAdjointDoubleLayerBoundaryOperator<
                        BasisFunctionType, ResultType>(
                            context, dualToRange, dualToRange /* hack */, domain, "ADLP"));
    }
};

template <typename BasisFunctionType, typename ResultType>
struct CreateOperator_K_plus_I_plus_5_V_plus_4_adj_Kp
{
    BoundaryOperator<BasisFunctionType, ResultType> operator() (
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange) const
    {
        return laplace3dDoubleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange, "DLP") +
                identityOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange, "Id") +
                5. * laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>(
                    context, domain, range, dualToRange, "SLP") +
                4. * adjoint(laplace3dAdjointDoubleLayerBoundaryOperator<
                        BasisFunctionType, ResultType>(
                            context, dualToRange, dualToRange /* hack */, domain, "ADLP"));
    }
};

template <typename BasisFunctionType, typename ResultType,
          template<typename BFT, typename RT> class CreateOperatorFunctor>
struct TestInDenseMode
{
    void test()
    {
        typedef ResultType RT;
        typedef typename ScalarTraits<ResultType>::RealType RealType;
        typedef BasisFunctionType BFT;

        CreateOperatorFunctor<BFT, RT> createOperator;

        GridParameters params;
        params.topology = GridParameters::TRIANGULAR;
        shared_ptr<Grid> grid = GridFactory::importGmshGrid(
                    params, "meshes/cube-12-reoriented.msh", false /* verbose */);

        shared_ptr<Space<BFT> > pwiseConstants(
                    new PiecewiseConstantScalarSpace<BFT>(grid));
        shared_ptr<Space<BFT> > pwiseLinears(
                    new PiecewiseLinearContinuousScalarSpace<BFT>(grid));

        AccuracyOptions accuracyOptions;
        accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
        shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                    new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

        AssemblyOptions assemblyOptionsRef;
        assemblyOptionsRef.setVerbosityLevel(VerbosityLevel::LOW);
        shared_ptr<Context<BFT, RT> > contextRef(
                    new Context<BFT, RT>(quadStrategy, assemblyOptionsRef));

        BoundaryOperator<BFT, RT> opRef = createOperator(
                    contextRef, pwiseLinears, pwiseLinears, pwiseConstants);
        Matrix<RT> weakFormRef = opRef.weakForm()->asMatrix();

        AssemblyOptions assemblyOptionsTest;
        assemblyOptionsTest.setVerbosityLevel(VerbosityLevel::LOW);
        assemblyOptionsTest.enableJointAssembly();
        shared_ptr<Context<BFT, RT> > contextTest(
                    new Context<BFT, RT>(quadStrategy, assemblyOptionsTest));

        BoundaryOperator<BFT, RT> opTest = createOperator(
                    contextTest, pwiseLinears, pwiseLinears, pwiseConstants);
        Matrix<RT> weakFormTest = opTest.weakForm()->asMatrix();

        BOOST_CHECK(check_arrays_are_close<RT>(
                        weakFormTest, weakFormRef,
                        10. * std::numeric_limits<RealType>::epsilon()));
    }
};

template <typename BasisFunctionType, typename ResultType,
          template<typename BFT, typename RT> class CreateOperatorFunctor>
struct TestInAcaMode
{
    void test()
    {
        typedef ResultType RT;
        typedef typename ScalarTraits<ResultType>::RealType RealType;
        typedef BasisFunctionType BFT;

        CreateOperatorFunctor<BFT, RT> createOperator;

        GridParameters params;
        params.topology = GridParameters::TRIANGULAR;
        shared_ptr<Grid> grid = GridFactory::importGmshGrid(
                    params, "meshes/cube-12-reoriented.msh", false /* verbose */);

        shared_ptr<Space<BFT> > pwiseConstants(
                    new PiecewiseConstantScalarSpace<BFT>(grid));
        shared_ptr<Space<BFT> > pwiseLinears(
                    new PiecewiseLinearContinuousScalarSpace<BFT>(grid));

        AccuracyOptions accuracyOptions;
        accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
        shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                    new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

        AcaOptions acaOptions;
        acaOptions.minimumBlockSize = 2;

        AssemblyOptions assemblyOptionsRef;
        assemblyOptionsRef.setVerbosityLevel(VerbosityLevel::LOW);
        assemblyOptionsRef.switchToAcaMode(acaOptions);
        assemblyOptionsRef.setMaxThreadCount(1);
        shared_ptr<Context<BFT, RT> > contextRef(
                    new Context<BFT, RT>(quadStrategy, assemblyOptionsRef));

        BoundaryOperator<BFT, RT> opRef = createOperator(
                    contextRef, pwiseLinears, pwiseLinears, pwiseConstants);
        Matrix<RT> weakFormRef = opRef.weakForm()->asMatrix();

        AssemblyOptions assemblyOptionsTest;
        assemblyOptionsTest.setVerbosityLevel(VerbosityLevel::LOW);
        assemblyOptionsTest.enableJointAssembly();
        assemblyOptionsTest.switchToAcaMode(acaOptions);
        assemblyOptionsTest.setMaxThreadCount(1);
        shared_ptr<Context<BFT, RT> > contextTest(
                    new Context<BFT, RT>(quadStrategy, assemblyOptionsTest));

        BoundaryOperator<BFT, RT> opTest = createOperator(
                    contextTest, pwiseLinears, pwiseLinears, pwiseConstants);
        Matrix<RT> weakFormTest = opTest.weakForm()->asMatrix();

        BOOST_CHECK(check_arrays_are_close<RT>(
                        weakFormTest, weakFormRef, 2. * acaOptions.eps));
    }
};

} // namespace

BOOST_AUTO_TEST_SUITE(JointAssembly)

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_V_plus_K,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_V_plus_K> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_3_V_plus_K,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_3_V_plus_K> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_V_plus_3_K,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_V_plus_3_K> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_V,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_V> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_3_V,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_3_V> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_V_plus_K_plus_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_V_plus_K_plus_Kp> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_3_V_plus_K_plus_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_3_V_plus_K_plus_Kp> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_V_plus_3_K_plus_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_V_plus_3_K_plus_Kp> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_I,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_I> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_3_I,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_3_I> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_V_plus_I,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_V_plus_I> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_3_I_plus_V,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_3_I_plus_V> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_I_plus_V,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_I_plus_V> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_I_plus_3_V,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_I_plus_3_V> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_2_I_plus_V_plus_I_plus_4_K,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_2_I_plus_V_plus_I_plus_4_K> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_K_plus_I_plus_4_adj_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_K_plus_I_plus_4_adj_Kp> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_K_plus_I_plus_5_V_plus_4_adj_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_K_plus_I_plus_5_V_plus_4_adj_Kp> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_dense_mode_for_K_plus_adj_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInDenseMode<BFT, RT, CreateOperator_K_plus_adj_Kp> test;
    test.test();
}

#ifdef WITH_AHMED

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_V_plus_K,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_V_plus_K> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_3_V_plus_K,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_3_V_plus_K> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_V_plus_3_K,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_V_plus_3_K> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_V,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_V> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_3_V,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_3_V> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_V_plus_K_plus_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_V_plus_K_plus_Kp> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_3_V_plus_K_plus_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_3_V_plus_K_plus_Kp> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_V_plus_3_K_plus_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_V_plus_3_K_plus_Kp> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_I,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_I> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_3_I,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_3_I> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_V_plus_I,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_V_plus_I> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_3_I_plus_V,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_3_I_plus_V> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_I_plus_V,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_I_plus_V> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_I_plus_3_V,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_I_plus_3_V> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_2_I_plus_V_plus_I_plus_4_K,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_2_I_plus_V_plus_I_plus_4_K> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_K_plus_adj_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_K_plus_adj_Kp> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_K_plus_I_plus_4_adj_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_K_plus_I_plus_4_adj_Kp> test;
    test.test();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(joint_assembly_works_in_aca_mode_for_K_plus_I_plus_5_V_plus_4_adj_Kp,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    TestInAcaMode<BFT, RT, CreateOperator_K_plus_I_plus_5_V_plus_4_adj_Kp> test;
    test.test();
}

#endif // WITH_AHMED

BOOST_AUTO_TEST_SUITE_END()

