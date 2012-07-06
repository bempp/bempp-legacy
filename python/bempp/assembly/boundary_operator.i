%{
#include "assembly/boundary_operator.hpp"
%}

// TODO
// %include "boundary_operator_docstrings.i"

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BoundaryOperator);

template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename BasisFunctionType> class Space;
class AssemblyOptions;
class Symmetry;

%extend BoundaryOperator
{
    %ignore BoundaryOperator;

    // SWIG doesn't handle properly the "using boost::shared_ptr" declaration
    // in BEM++'s code, therefore these wrappers returning fully qualified
    // shared_ptrs are needed
    boost::shared_ptr<const DiscreteBoundaryOperator<ResultType> > weakForm() const
    {
        return $self->weakForm();
    }
    %ignore weakForm;

    boost::shared_ptr<const AbstractBoundaryOperator<BasisFunctionType, ResultType> > abstractOperator() const
    {
        return $self->abstractOperator();
    }
    %ignore abstractOperator;

    boost::shared_ptr<const Context<BasisFunctionType, ResultType> > context() const
    {
        return $self->context();
    }
    %ignore context;

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

    BoundaryOperator<BasisFunctionType, ResultType> __add__(
        const BoundaryOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self + other;
    }

    BoundaryOperator<BasisFunctionType, ResultType> __sub__(
        const BoundaryOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self - other;
    }

    BoundaryOperator<BasisFunctionType, ResultType> __mul__(
        ResultType other)
    {
        return *$self * other;
    }

    GridFunction<BasisFunctionType, ResultType> __mul__(
        const GridFunction<BasisFunctionType, ResultType>& other)
    {
        return *$self * other;
    }

    BoundaryOperator<BasisFunctionType, ResultType> __rmul__(
        ResultType other)
    {
        return *$self * other;
    }

    BoundaryOperator<BasisFunctionType, ResultType> __div__(
        ResultType other)
    {
        return *$self / other;
    }
}

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BoundaryOperator);

} // namespace Bempp

%include "assembly/boundary_operator.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(BoundaryOperator);
}
