#ifndef bempp_abstract_boundary_operator_id_hpp
#define bempp_abstract_boundary_operator_id_hpp

namespace Bempp
{

class AbstractBoundaryOperatorId
{
public:
    virtual ~AbstractBoundaryOperatorId() {}
    virtual bool less(const AbstractBoundaryOperatorId& other) const = 0;
};

} // namespace Bempp

#endif // bempp_abstract_boundary_operator_id_hpp
