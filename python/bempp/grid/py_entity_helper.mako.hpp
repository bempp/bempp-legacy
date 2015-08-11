#ifndef bempp_py_entity_helper_hpp
#define bempp_py_entity_helper_hpp


#include "bempp/grid/entity.hpp"
#include "bempp/grid/entity_pointer.hpp"

namespace Bempp {

    inline std::unique_ptr<EntityPointer<0>> py_get_father_from_entity(const Entity<0>& entity)
    {
        return entity.father();
    }

}

#endif
