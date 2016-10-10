#ifndef bempp_p1_entity_iterator_hpp
#define bempp_p1_entity_iterator_hpp

#include "p1_entity_pointer.hpp"

namespace BemppGrid {

template <int codim>
class P1EntityIteratorImp : public EntityPointerImp<codim>
{

    public:

        virtual void increment();

};


}


#endif
