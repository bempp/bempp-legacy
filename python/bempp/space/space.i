%{
  #include "space/space.hpp"
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

namespace Bempp{
  %template(SpaceDouble) Space<double>;
     }
