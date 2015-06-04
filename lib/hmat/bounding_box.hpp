#ifndef HMAT_BOUNDING_BOX_HPP
#define HMAT_BOUNDING_BOX_HPP

#include "common.hpp"
#include <array>
#include <vector>
#include <iostream>

namespace hmat {

class BoundingBox {

public:
  BoundingBox();

  BoundingBox(double xmin_, double xmax_, double ymin_, double ymax_,
              double zmin_, double zmax_);

  BoundingBox(const std::array<double, 6> &bounds);

  BoundingBox(const BoundingBox &other);
  BoundingBox(const BoundingBox &&other);

  BoundingBox &operator=(const BoundingBox &other);
  BoundingBox &operator=(const BoundingBox &&other);

  double xmin() const;
  double xmax() const;
  double ymin() const;
  double ymax() const;
  double zmin() const;
  double zmax() const;

  std::pair<BoundingBox, BoundingBox> divide(int dim, double f) const;
  std::pair<BoundingBox, BoundingBox> divideMaxDimension(double f) const;
  void merge(const BoundingBox &other);
  bool contains(const std::array<double, 3> point);
  const std::array<double, 6> &bounds() const;

  const std::array<double, 3> cornerPoint(int index) const;

  void corners(std::vector<Point> &points) const;

  int maxDimension() const;
  double diameter() const;

  double distance(const BoundingBox &other) const;

private:
  void computeGeometricData();
  double intervalDistance(double a1, double b1, double a2, double b2) const;

  std::array<double, 6> m_boundingBoxData;
  int m_maxDimension;
  double m_diameter;
};

BoundingBox boundingBoxFromTriangleData(std::array<double, 3> p1,
                                        std::array<double, 3> p2,
                                        std::array<double, 3> p3);

std::ostream &operator<<(std::ostream &os, const BoundingBox &box);
}

#include "bounding_box_impl.hpp"

#endif
