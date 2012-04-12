#ifndef fiber_quadrature_options_hpp
#define fiber_quadrature_options_hpp

namespace Fiber
{

struct QuadratureOptions
{
    QuadratureOptions() : mode(ORDER_INCREMENT), orderIncrement(0) {
    }

    enum Mode {
        ORDER_INCREMENT, EXACT_ORDER
    } mode;

    union {
        int orderIncrement;
        int order;
    };
};

}

#endif
