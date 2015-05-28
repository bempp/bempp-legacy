// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_BOUNDING_BOX_IMPL_HPP
#define HMAT_BOUNDING_BOX_IMPL_HPP

#include "bounding_box.hpp"
#include <cmath>
#include <cassert>

namespace hmat {

inline BoundingBox::BoundingBox()
    : BoundingBox(std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::min(),
                  std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::min(),
                  std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::min()) {}

inline BoundingBox::BoundingBox(double xmin_, double xmax_, double ymin_,
                                double ymax_, double zmin_, double zmax_)
    : BoundingBox(
          std::array<double, 6>{{xmin_, xmax_, ymin_, ymax_, zmin_, zmax_}}) {}

inline BoundingBox::BoundingBox(const std::array<double, 6> &bounds)
    : m_boundingBoxData(bounds) {
  computeGeometricData();
}

inline BoundingBox::BoundingBox(const BoundingBox &other)
    : m_boundingBoxData(other.m_boundingBoxData),
      m_maxDimension(other.m_maxDimension), m_diameter(other.m_diameter) {}

inline BoundingBox::BoundingBox(const BoundingBox &&other)
    : m_boundingBoxData(std::move(other.m_boundingBoxData)),
      m_maxDimension(other.m_maxDimension), m_diameter(other.m_diameter) {}

inline BoundingBox &BoundingBox::operator=(const BoundingBox &other) {
  m_boundingBoxData = other.m_boundingBoxData;
  m_maxDimension = other.m_maxDimension;
  m_diameter = other.m_diameter;
  return *this;
}

inline BoundingBox &BoundingBox::operator=(const BoundingBox &&other) {
  m_boundingBoxData = std::move(other.m_boundingBoxData);
  m_maxDimension = other.m_maxDimension;
  m_diameter = other.m_diameter;
  return *this;
}

inline void BoundingBox::computeGeometricData() {

  std::array<double, 3> diameters;
  diameters[0] = m_boundingBoxData[1] - m_boundingBoxData[0];
  diameters[1] = m_boundingBoxData[3] - m_boundingBoxData[2];
  diameters[2] = m_boundingBoxData[5] - m_boundingBoxData[4];

  m_maxDimension = 0;
  m_diameter = 0;

  for (int i = 0; i < 3; ++i) {
    if (diameters[i] > diameters[m_maxDimension]) {
      m_maxDimension = i;
    }
  }
  m_diameter =
      std::sqrt(diameters[0] * diameters[0] + diameters[1] * diameters[1] +
                diameters[2] * diameters[2]);
}

inline double BoundingBox::xmin() const { return m_boundingBoxData[0]; }

inline double BoundingBox::xmax() const { return m_boundingBoxData[1]; }

inline double BoundingBox::ymin() const { return m_boundingBoxData[2]; }

inline double BoundingBox::ymax() const { return m_boundingBoxData[3]; }

inline double BoundingBox::zmin() const { return m_boundingBoxData[4]; }

inline double BoundingBox::zmax() const { return m_boundingBoxData[5]; }

inline void BoundingBox::merge(const BoundingBox &other) {

  m_boundingBoxData[0] =
      std::min(m_boundingBoxData[0], other.m_boundingBoxData[0]);
  m_boundingBoxData[1] =
      std::max(m_boundingBoxData[1], other.m_boundingBoxData[1]);

  m_boundingBoxData[2] =
      std::min(m_boundingBoxData[2], other.m_boundingBoxData[2]);
  m_boundingBoxData[3] =
      std::max(m_boundingBoxData[3], other.m_boundingBoxData[3]);

  m_boundingBoxData[4] =
      std::min(m_boundingBoxData[4], other.m_boundingBoxData[4]);
  m_boundingBoxData[5] =
      std::max(m_boundingBoxData[5], other.m_boundingBoxData[5]);

  computeGeometricData();
}

inline std::pair<BoundingBox, BoundingBox> BoundingBox::divide(int dim,
                                                               double f) const {
  std::pair<BoundingBox, BoundingBox> result{*this, *this};

  assert(0 <= dim && dim <= 2);

  if (dim == 0) {
    result.first.m_boundingBoxData[1] =
        m_boundingBoxData[0] +
        f * (m_boundingBoxData[1] - m_boundingBoxData[0]);
    result.second.m_boundingBoxData[0] = result.first.m_boundingBoxData[1];

  } else if (dim == 1) {
    result.first.m_boundingBoxData[3] =
        m_boundingBoxData[2] +
        f * (m_boundingBoxData[3] - m_boundingBoxData[2]);
    result.second.m_boundingBoxData[2] = result.first.m_boundingBoxData[3];

  } else if (dim == 2) {
    result.first.m_boundingBoxData[5] =
        m_boundingBoxData[4] +
        f * (m_boundingBoxData[5] - m_boundingBoxData[4]);
    result.second.m_boundingBoxData[4] = result.first.m_boundingBoxData[5];

  } else
    throw std::invalid_argument("BoundingBox::divide: wrong index.");

  result.first.computeGeometricData();
  result.second.computeGeometricData();

  return result;
}

inline std::pair<BoundingBox, BoundingBox>
BoundingBox::divideMaxDimension(double f) const {
  return divide(m_maxDimension, f);
}

inline bool BoundingBox::contains(const std::array<double, 3> point) {

  if (m_boundingBoxData[0] <= point[0] && m_boundingBoxData[1] > point[0] &&
      m_boundingBoxData[2] <= point[1] && m_boundingBoxData[3] > point[1] &&
      m_boundingBoxData[4] <= point[2] && m_boundingBoxData[5] > point[2])
    return true;
  else
    return false;
}

inline const std::array<double, 6> &BoundingBox::bounds() const {

  return m_boundingBoxData;
}

inline int BoundingBox::maxDimension() const { return m_maxDimension; }

inline double BoundingBox::diameter() const { return m_diameter; }

inline const std::array<double, 3> BoundingBox::cornerPoint(int index) const {

  assert(0 <= index && index <= 7);

  std::array<double, 3> result = {{0, 0, 0}};

  switch (index) {
  case 0:
    result[0] = xmin();
    result[1] = ymin();
    result[2] = zmin();
    break;
  case 1:
    result[0] = xmax();
    result[1] = ymin();
    result[2] = zmin();
    break;
  case 2:
    result[0] = xmax();
    result[1] = ymax();
    result[2] = zmin();
    break;
  case 3:
    result[0] = xmin();
    result[1] = ymax();
    result[2] = zmin();
    break;
  case 4:
    result[0] = xmin();
    result[1] = ymin();
    result[2] = zmax();
    break;
  case 5:
    result[0] = xmax();
    result[1] = ymin();
    result[2] = zmax();
    break;
  case 6:
    result[0] = xmax();
    result[1] = ymax();
    result[2] = zmax();
    break;
  case 7:
    result[0] = xmin();
    result[1] = ymax();
    result[2] = zmax();
    break;
  }
  return result;
}
  
inline const std::vector<Point> 
BoundingBox::corners() const {

  std::vector<Point> result;
  result.reserve(8);
  for (int i = 0; i < 8; ++i)
  {
    std::array<double,3> corner = cornerPoint(i);
    result.push_back(Point(corner[0],corner[1],corner[2]));
  }
  return result;

}

inline double BoundingBox::distance(const BoundingBox &other) const {
  auto bounds = other.bounds();

  double sum = 0;
  for (int i = 0; i < 3; ++i) {
    double dist =
        intervalDistance(m_boundingBoxData[2 * i], m_boundingBoxData[2 * i + 1],
                         bounds[2 * i], bounds[2 * i + 1]);
    sum += dist * dist;
  }

  return std::sqrt(sum);
}

inline double BoundingBox::intervalDistance(double a1, double b1, double a2,
                                            double b2) const {
  if (a2 > b1)
    return a2 - b1;
  else if (a1 > b2)
    return a1 - b2;
  else
    return 0;
}

inline std::ostream &operator<<(std::ostream &os, const BoundingBox &box) {
  auto bounds = box.bounds();

  for (auto it = begin(bounds); it != end(bounds); ++it)
    os << *it << " ";

  return os;
}

inline BoundingBox boundingBoxFromTriangleData(std::array<double, 3> p1,
                                               std::array<double, 3> p2,
                                               std::array<double, 3> p3) {

  double xmin = std::min({p1[0], p2[0], p3[0]});
  double xmax = std::max({p1[0], p2[0], p3[0]});

  double ymin = std::min({p1[1], p2[1], p3[1]});
  double ymax = std::max({p1[1], p2[1], p3[1]});

  double zmin = std::min({p1[2], p2[2], p3[2]});
  double zmax = std::max({p1[2], p2[2], p3[2]});

  return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
}
}

#endif
