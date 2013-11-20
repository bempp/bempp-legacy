#ifndef bempp_type_combinations_hpp
#define bempp_type_combinations_hpp

#include "bempp/common/config_data_types.hpp"

#include <boost/mpl/list.hpp>
#include <complex>

template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
struct BasisKernelResultTypeTraits
{
    typedef BasisFunctionType_ BasisFunctionType;
    typedef KernelType_ KernelType;
    typedef ResultType_ ResultType;
};

#if defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#  if defined(ENABLE_COMPLEX_KERNELS)

typedef boost::mpl::list<
#    if defined(ENABLE_SINGLE_PRECISION)
BasisKernelResultTypeTraits<float, float, float>,
BasisKernelResultTypeTraits<float, float, std::complex<float> >,
BasisKernelResultTypeTraits<float, std::complex<float>, std::complex<float> >,
BasisKernelResultTypeTraits<std::complex<float>, float, std::complex<float> >,
BasisKernelResultTypeTraits<std::complex<float>, std::complex<float>, std::complex<float> >
#    endif
#    if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_DOUBLE_PRECISION)
,
#    endif
#    if defined(ENABLE_DOUBLE_PRECISION)
BasisKernelResultTypeTraits<double, double, double>,
BasisKernelResultTypeTraits<double, double, std::complex<double> >,
BasisKernelResultTypeTraits<double, std::complex<double>, std::complex<double> >,
BasisKernelResultTypeTraits<std::complex<double>, double, std::complex<double> >,
BasisKernelResultTypeTraits<std::complex<double>, std::complex<double>, std::complex<double> >
#    endif
>
basis_kernel_result_combinations;

#  else // !defined(ENABLE_COMPLEX_KERNELS)

typedef boost::mpl::list<
#    if defined(ENABLE_SINGLE_PRECISION)
BasisKernelResultTypeTraits<float, float, float>,
BasisKernelResultTypeTraits<float, float, std::complex<float> >,
BasisKernelResultTypeTraits<std::complex<float>, float, std::complex<float> >
#    endif
#    if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_DOUBLE_PRECISION)
,
#    endif
#    if defined(ENABLE_DOUBLE_PRECISION)
BasisKernelResultTypeTraits<double, double, double>,
BasisKernelResultTypeTraits<double, double, std::complex<double> >,
BasisKernelResultTypeTraits<std::complex<double>, double, std::complex<double> >
#    endif
>
basis_kernel_result_combinations;

#  endif
#else // !defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#  if defined(ENABLE_COMPLEX_KERNELS)

typedef boost::mpl::list<
#    if defined(ENABLE_SINGLE_PRECISION)
BasisKernelResultTypeTraits<float, float, float>,
BasisKernelResultTypeTraits<float, float, std::complex<float> >,
BasisKernelResultTypeTraits<float, std::complex<float>, std::complex<float> >
#    endif
#    if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_DOUBLE_PRECISION)
,
#    endif
#    if defined(ENABLE_DOUBLE_PRECISION)
BasisKernelResultTypeTraits<double, double, double>,
BasisKernelResultTypeTraits<double, double, std::complex<double> >,
BasisKernelResultTypeTraits<double, std::complex<double>, std::complex<double> >
#    endif
>
basis_kernel_result_combinations;

#  else // !defined(ENABLE_COMPLEX_KERNELS)

typedef boost::mpl::list<
#    if defined(ENABLE_SINGLE_PRECISION)
BasisKernelResultTypeTraits<float, float, float>
#    endif
#    if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_DOUBLE_PRECISION)
,
#    endif
#    if defined(ENABLE_DOUBLE_PRECISION)
BasisKernelResultTypeTraits<double, double, double>
#    endif
>
basis_kernel_result_combinations;

#  endif
#endif


template <typename BasisFunctionType_, typename ResultType_>
struct BasisResultTypeTraits
{
    typedef BasisFunctionType_ BasisFunctionType;
    typedef ResultType_ ResultType;
};

#if defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)

typedef boost::mpl::list<
#  if defined(ENABLE_SINGLE_PRECISION)
BasisResultTypeTraits<float, float>,
BasisResultTypeTraits<float, std::complex<float> >,
BasisResultTypeTraits<std::complex<float>, std::complex<float> >
#  endif
#  if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_DOUBLE_PRECISION)
,
#  endif
#  if defined(ENABLE_DOUBLE_PRECISION)
BasisResultTypeTraits<double, double>,
BasisResultTypeTraits<double, std::complex<double> >,
BasisResultTypeTraits<std::complex<double>, std::complex<double> >
#  endif
>
basis_result_combinations;

#else // !defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#  if defined(ENABLE_COMPLEX_KERNELS)

typedef boost::mpl::list<
#    if defined(ENABLE_SINGLE_PRECISION)
BasisResultTypeTraits<float, float>,
BasisResultTypeTraits<float, std::complex<float> >
#    endif
#    if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_DOUBLE_PRECISION)
,
#    endif
#    if defined(ENABLE_DOUBLE_PRECISION)
BasisResultTypeTraits<double, double>,
BasisResultTypeTraits<double, std::complex<double> >
#    endif
>
basis_result_combinations;

#  else // !defined(ENABLE_COMPLEX_KERNELS)

typedef boost::mpl::list<
#    if defined(ENABLE_SINGLE_PRECISION)
BasisResultTypeTraits<float, float>
#    endif
#    if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_DOUBLE_PRECISION)
,
#    endif
#    if defined(ENABLE_DOUBLE_PRECISION)
BasisResultTypeTraits<double, double>
#    endif
>
basis_result_combinations;

#  endif
#endif

#endif
