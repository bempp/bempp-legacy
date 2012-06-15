%{
  #include "space/space.hpp"
  #include <complex>
%}

// TODO
// %include "space_docstrings.i"



namespace Bempp 
{

  template<typename BasisFunctionType> class Space;


  %extend Space
  {

    // this function is only for internal use
    %ignore basis;

    // this function is only for internal use
    %ignore shapeFunctionValueExpression;

    // to be wrapped later...
    %ignore setElementVariant;

    // to be wrapped later...
    %ignore elementVariant;


    %ignore globalDofs;

    %ignore global2localDofs;

    // this function is only for internal use
    %ignore globalDofPositions;


    %ignore dumpClusterIds;

    %ignore applyMassMatrix;

    %ignore applyInverseMassMatrix;

  }



} // namespace Bempp

%include "space/space.hpp"

namespace Bempp
{

  %template(Space_float64)     Space<double>;
  %template(Space_float32)     Space<float>;
  %template(Space_complex64)   Space<std::complex<float> >;
  %template(Space_complex128)  Space<std::complex<double> >;

}

