%{
#include "assembly/abstract_boundary_operator.hpp"
#include "assembly/abstract_boundary_operator_sum.hpp"
#include <complex>
%}

// TODO
// %include "abstract_boundary_operator_docstrings.i"

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);

template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename BasisFunctionType> class Space;
class AssemblyOptions;
class Symmetry;

%extend AbstractBoundaryOperator
{
    // SWIG doesn't handle properly the "using boost::shared_ptr" declaration
    // in BEM++'s code, therefore these wrappers returning fully qualified
    // shared_ptrs are needed
    boost::shared_ptr<const DiscreteBoundaryOperator<ResultType> >
    assembleWeakForm(const Context<BasisFunctionType, ResultType>& context) const
    {
        return $self->assembleWeakForm(context);
    }
    %ignore assembleWeakForm;

    boost::shared_ptr<const Space<BasisFunctionType> > domain() const
    {
        return $self->domain();
    }
    %ignore domain;

    boost::shared_ptr<const Space<BasisFunctionType> > range() const
    {
        return $self->range();
    }
    %ignore range;

    boost::shared_ptr<const Space<BasisFunctionType> > dualToRange() const
    {
        return $self->dualToRange();
    }
    %ignore dualToRange;

    %ignore id;
}

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);

} // namespace Bempp

%include "assembly/abstract_boundary_operator.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);
}
