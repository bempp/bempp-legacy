#ifndef HMAT_POINT_HPP
#define HMAT_POINT_HPP

#include "common.hpp"

namespace hmat {

class Point
{
public:
    Point();
    Point(double x, double y, double z);

    const double& x() const;
    double& x();

    const double& y() const;
    double& y();

    const double& z() const;
    double& z();

private:
    double m_x;
    double m_y;
    double m_z;

};
}

#include "point_impl.hpp"
#endif
