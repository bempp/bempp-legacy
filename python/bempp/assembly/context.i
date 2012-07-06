%{
#include "assembly/context.hpp"
%}

// TODO
// %include "context_docstrings.i"

%shared_ptr(Bempp::Context<float, float>);
%shared_ptr(Bempp::Context<float, std::complex<float> >);
%shared_ptr(Bempp::Context<std::complex<float>, std::complex<float> >);
%shared_ptr(Bempp::Context<double, double>);
%shared_ptr(Bempp::Context<double, std::complex<double> >);
%shared_ptr(Bempp::Context<std::complex<double>, std::complex<double> >);

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Context);

class AssemblyOptions;

%extend Context
{

    Context(const boost::shared_ptr<typename Context<BasisFunctionType, 
            ResultType>::LocalAssemblerFactory>& localAssemblerFactory,
            const AssemblyOptions& assemblyOptions) 
    {
        return new Bempp::Context<BasisFunctionType, ResultType >(localAssemblerFactory, assemblyOptions);
    }
    %ignore Context;

    boost::shared_ptr<const DiscreteBoundaryOperator<ResultType> >
    getWeakForm(const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op) const {
        return $self->getWeakForm(op);
    }
    %ignore getWeakForm;

}

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Context);

} // namespace Bempp

%include "assembly/context.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Context);
}
