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

    // this function is only for internal use
    %ignore setElementVariant;

    // this function is only for internal use
    %ignore elementVariant;


    %ignore globalDofs;

    %ignore global2localDofs;

    // this function is only for internal use
    %ignore globalDofPositions;
    
    // this function is only for internal use
    %ignore globalDofPositions;


    %ignore applyMassMatrix;
    %ignore applyInverseMassMatrix;

    %ignore dumpClusterIds;

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

%pythoncode %{

class Space(object):

  def assignDofs(self):
      """Assign degrees of freedom to the space"""
      pass

     %}
