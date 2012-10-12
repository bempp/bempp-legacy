#ifndef hermite_interpolator_hpp
#define hermite_interpolator_hpp

#include "../common/common.hpp"
#include "scalar_traits.hpp"

#include <cassert>
#include <stdexcept>
#include <vector>

namespace Fiber
{

template <typename ValueType>
class HermiteInterpolator
{
public:
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    HermiteInterpolator() : m_start(0.), m_end(0.), m_n(0), m_interval(0.)
        {}

    CoordinateType rangeStart() { return m_start; }
    CoordinateType rangeEnd() { return m_end; }

    void initialize(CoordinateType start, CoordinateType end,
                    const std::vector<ValueType>& values,
                    const std::vector<ValueType>& derivatives) {
        if (values.size() != derivatives.size())
            throw std::invalid_argument("HermiteInterpolator::setData(): "
                                        "'values' and 'derivatives' must "
                                        "have the same length");
        if (values.size() < 2)
            throw std::invalid_argument("HermiteInterpolator::setData(): "
                                        "at least two points are required");
        if (end <= start)
            throw std::invalid_argument("HermiteInterpolator::setData(): "
                                        "'start' must be smaller than 'end'");
        m_start = start;
        m_end = end;
        m_n = values.size();
        m_interval = (end - start) / (m_n - 1);
        m_values = values;
        m_derivatives = derivatives;
    }

    ValueType evaluate(CoordinateType x) const {
        assert(x >= m_start && x <= m_end);
        CoordinateType fpn;
        const CoordinateType t = std::modf((x - m_start) / m_interval, &fpn);
        const int n = int(fpn);
        assert(n >= 0 && n < m_n);
        assert(t >= 0 && t <= 1);
        const CoordinateType one = 1., two = 2., three = 3.;
        const CoordinateType tm1 = t - one;
        // std::cout << x << n << " " << t << " " << m_interval << " " << m_values[n] << " " << m_derivatives[n] << " " << "\n";
        return (m_values[n] * (one + two * t) * tm1 * tm1 +
                m_interval * m_derivatives[n] * t * tm1 * tm1 +
                m_values[n+1] * t * t * (three - two * t) +
                m_interval * m_derivatives[n+1] * t * t * tm1);
    }

private:
    /** \cond PRIVATE */
    CoordinateType m_start, m_end;
    int m_n;
    CoordinateType m_interval;
    std::vector<ValueType> m_values;
    std::vector<ValueType> m_derivatives;
    /** \endcond */
};

} // namespace Fiber

#endif
