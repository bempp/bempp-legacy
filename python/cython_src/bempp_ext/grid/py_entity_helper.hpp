#ifndef bempp_py_entity_helper_hpp
#define bempp_py_entity_helper_hpp


#include "bempp/grid/entity.hpp"
#include "bempp/grid/entity_pointer.hpp"
#include "bempp/grid/entity_iterator.hpp"

namespace Bempp {

    inline std::unique_ptr<EntityPointer<0>> py_get_father_from_entity(const Entity<0>& entity)
    {
        return entity.father();
    }

    template <int iteratorCodim, int entityCodim>
    inline std::unique_ptr<EntityIterator<iteratorCodim>>
    py_get_entity_iterator_from_entity(const Entity<entityCodim>& entity)
    {
        //return std::unique_ptr<EntityIterator<iteratorCodim>>(nullptr);
        return entity.template subEntityIterator<iteratorCodim>();

    }

}

#endif
