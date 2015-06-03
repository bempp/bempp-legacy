// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_COMMON_HPP
#define HMAT_COMMON_HPP

#include "shared_ptr.hpp"

#include <CGAL/Simple_cartesian.h>

#include <array>
#include <vector>

namespace hmat {

typedef std::array<std::size_t, 4> BlockIndexRangeType;
typedef std::array<std::size_t, 2> IndexRangeType;
typedef std::vector<std::size_t> IndexSetType;

typedef CGAL::Simple_cartesian<double> CgalKernel;
typedef CgalKernel::Point_3 Point;
typedef CgalKernel::Line_3 Line;
typedef CgalKernel::Plane_3 Plane;

enum RowColSelector { ROW, COL };

enum TransposeMode { NOTRANS, TRANS, CONJ, CONJTRANS };

enum DataBlockType {

  DENSE,
  LOW_RANK_AB

};

IndexSetType fillIndexRange(std::size_t start, std::size_t stop);

double obbDistance(const std::vector<Point> &points1,
                   const std::vector<Point> &points2);

double clusterDiameter(const std::vector<Point> &points);

std::vector<Point> convex_hull(const std::vector<Point> &points);
}

#include "common_impl.hpp"

#endif
