%{
#include "linalg/preconditioner.hpp"
%}



namespace std {

%template(vector_DiscreteOperator_float32) vector<boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<float> > >;
%template(vector_DiscreteOperator_float64) vector<boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<double> > >;
%template(vector_DiscreteOperatorVector_complex64) vector<boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<std::complex<float> > > >;
%template(vector_DiscreteOperatorVector_complex128) vector<boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<std::complex<double> > > >;



}



namespace Bempp {



BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(Preconditioner)
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(acaDiscreteOperatorToPreconditioner)
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(acaBlockDiagonalPreconditioner)
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(discreteBlockDiagonalPreconditioner)


%extend Preconditioner {

    %ignore Preconditioner;
    %ignore get;
    %ignore trilinosPreconditioner;
}

BEMPP_EXTEND_CLASS_TEMPLATED_ON_VALUE(Preconditioner)
}

#define shared_ptr boost::shared_ptr
%include "linalg/preconditioner.hpp"
#undef shared_ptr

namespace Bempp {
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(Preconditioner)
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(acaDiscreteOperatorToPreconditioner)
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(acaBlockDiagonalPreconditioner)
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(discreteBlockDiagonalPreconditioner)

}

