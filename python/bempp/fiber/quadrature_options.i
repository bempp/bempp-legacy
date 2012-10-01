%{
#include "fiber/quadrature_options.hpp"
%}

%extend Fiber::QuadratureOptions
{
    %pythoncode %{
    def _setOrderIncrement(self, value):
        print ("Warning: the orderIncrement property is deprecated and "
               "maintained only for backwards compatibility. Please use "
               "setRelativeQuadratureOrder(quadratureOrder) instead.")
        self.setRelativeQuadratureOrder(value)

    def _setOrder(self, value):
        print ("Warning: the order property is deprecated and "
               "maintained only for backwards compatibility. Please use "
               "setAbsoluteQuadratureOrder(quadratureOrder) instead.")
        self.setAbsoluteQuadratureOrder(value)

    def _orderOrOrderIncrement(self):
        return self.quadratureOrder(0)

    orderIncrement = property(_orderOrOrderIncrement, _setOrderIncrement,
        """Relative quadrature order.

           DEPRECATED. Please use setRelativeQuadratureOrder() instead.
        """)
    order = property(_orderOrOrderIncrement, _setOrder,
        """Relative quadrature order.

        DEPRECATED. Please use setAbsoluteQuadratureOrder() instead.
        """)
    %}
}

%include "fiber/quadrature_options.hpp"
