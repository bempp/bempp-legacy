%{
#include "linalg/solver.hpp"
%}

%include "std_vector.i"
namespace std
{
%template(vector_GridFunction_float32_float32)
    vector<Bempp::GridFunction<float, float> >;
%template(vector_GridFunction_float32_complex64)
    vector<Bempp::GridFunction<float, std::complex<float> > >;
%template(vector_GridFunction_complex64_complex64)
    vector<Bempp::GridFunction<std::complex<float>, std::complex<float> > >;

%template(vector_GridFunction_float64_float64)
    vector<Bempp::GridFunction<double, double> >;
%template(vector_GridFunction_float64_complex128)
    vector<Bempp::GridFunction<double, std::complex<double> > >;
%template(vector_GridFunction_complex128_complex128)
    vector<Bempp::GridFunction<std::complex<double>, std::complex<double> > >;
}

namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Solver);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Solver);
}

#define shared_ptr boost::shared_ptr
%include "linalg/solver.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Solver);
}
