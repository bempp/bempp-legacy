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

%pythoncode %{

class Space(object):

  def domainDimension(self):
      """Dimension of the surface on which the functions are defined."""
      pass

  def codomainDimension(self):
      """Dimension of the codomain of the functions.

         In other words, number of components of the values of the functions.
         (E.g. H1 space -> 1, H(curl) space on a 2D surface -> 2).
         virtual int codomainDimension() const = 0;
      """
      pass

  def grid(self):
      """Reference to the grid on which the functions are defined."""
      pass

  def assignDofs(self):
      """Assign global degrees of freedom to local degrees of freedom."""
      pass

  def dofsAssigned(self):
      """True if assignDofs() has been called before, false otherwise."""
      pass

  def globalDofCount(self):
      """Number of global degrees of freedom."""
      pass     

  def global2localDofs(self):
      pass

     %}
