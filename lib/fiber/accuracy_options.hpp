#ifndef fiber_accuracy_options_hpp
#define fiber_accuracy_options_hpp

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

class AccuracyOptions
{
public:
    QuadratureOptions regular;
    QuadratureOptions singular;
};

}

#endif
