#ifndef BEMPP_PYTHON_SPACE_INPLACE_H
#define BEMPP_PYTHON_SPACE_INPLACE_H

#include "bempp/assembly/boundary_operator.hpp"

namespace {

    // Constructs C++ object in python-managed memory
    template<class BASIS, class RESULT>
        inline void inplace_boundary_operator(void *_memory) {
            new(_memory) Bempp::BoundaryOperator<BASIS, RESULT>();
        }

    // Deconstructs C++ object
    template<class BASIS, class RESULT>
        inline void deconstructor(void *_memory) {
            static_cast< Bempp::BoundaryOperator<BASIS, RESULT>* >(_memory)
                ->~BoundaryOperator<BASIS, RESULT>();
        }
}
#endif
