#ifndef HMAT_POINT_IMPL_HPP
#define HMAT_POINT_IMPL_HPP


#include "point.hpp"

namespace hmat {

inline Point::Point():
    m_x(0), m_y(0), m_z(0) {}

inline Point::Point(double x, double y, double z):
    m_x(x), m_y(y), m_z(z){}

inline const double& Point::x() const
{
    return m_x;

}

inline const double& Point::y() const
{
    return m_y;

}

inline const double& Point::z() const
{
    return m_z;
}

inline double& Point::x() 
{
    return m_x;
}

inline double& Point::y()
{
    return m_y;
}

inline double& Point::z()
{
    return m_z;
}

}

#endif
