#ifndef BEMPP_PYTHON_SPACE_INPLACE_H
#define BEMPP_PYTHON_SPACE_INPLACE_H

#include "bempp/assembly/boundary_operator.hpp"

namespace {

    template<class BASIS, class RESULT>
        inline void inplace_boundary_operator(void *_memory) {
            new(_memory) Bempp::BoundaryOperator<BASIS, RESULT>();
        }
}
#endif
