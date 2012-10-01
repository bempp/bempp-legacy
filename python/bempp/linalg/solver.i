%{
#include "linalg/solver.hpp"
%}


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

namespace Bempp {

%typemap(in) ConvergenceTestMode::Mode
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError,
                        "in method '$symname', argument $argnum: "
                        "expected a string");
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "test_convergence_in_range")
        $1 = Bempp::ConvergenceTestMode::TEST_CONVERGENCE_IN_RANGE;
    else if (s == "test_convergence_in_dual_to_range")
        $1 = Bempp::ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE;
    else
    {
        PyErr_SetString(PyExc_ValueError,
                        "in method '$symname', argument $argnum: "
                        "expected one of 'test_convergence_in_range' or "
                        "'test_convergence_in_dual_to_range'.");
        SWIG_fail;
    }
}

%typemap(typecheck) ConvergenceTestMode::Mode = std::string;

}


namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Solver);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Solver);

// Convert any empty elements in the supplied list to uninitialized
// GridFunctions
%define BEMPP_EXTEND_SOLVER(BASIS, RESULT, PYBASIS, PYRESULT)
    %extend Solver< BASIS, RESULT >
    {
        %pythonprepend solve %{
            try:
                newargs = ([],)
                for a in args[0]:
                    if a:
                        newargs[0].append(a)
                    else:
                        newargs[0].append(globals()[
                            "uninitializedGridFunction_"+
                            self.basisFunctionType()+"_"+self.resultType()]())
                args = newargs
            except: # keep the list of arguments unchanged
                pass
        %}
    }
%enddef
BEMPP_ITERATE_OVER_BASIS_AND_RESULT_TYPES(BEMPP_EXTEND_SOLVER)

}

#define shared_ptr boost::shared_ptr
%include "linalg/solver.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Solver);
}
