/*
 *                This source code is part of
 * 
 *                     E  R  K  A  L  E
 *                             -
 *                       DFT from Hel
 *
 * Written by Susi Lehtola, 2010-2011
 * Copyright (c) 2010-2011, Susi Lehtola
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 */



#ifndef ERKALE_LEBEDEV
#define ERKALE_LEBEDEV

#include <vector>

/// Orders of supported rules (integrates spherical harmonics exactly up to given order)
const int lebedev_orders[]={3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131};/// The number of points used by the supported rules
const int lebedev_degrees[]={6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266,302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};

/// Coordinates and weight for a point in Lebedev quadrature
typedef struct {
  /// x coordinate of point
  double x;
  /// y coordinate of point
  double y;
  /// z coordinate of point
  double z;
  /// Angular weight of point
  double w;
} lebedev_point_t;

/// Get a Lebedev sphere of the wanted order.
std::vector<lebedev_point_t> lebedev_sphere(int order);
/// Determine the next order of quadrature, if one is supported.
int next_lebedev(int order);

/// Worker routine - get a Lebedev sphere with the wanted amount of points
std::vector<lebedev_point_t> getLebedevSphere(int npoints);
/// Worker routine - add points to the sphere structure
void getLebedevReccurencePoints(int type, int & start, double a, double b, double v, std::vector<lebedev_point_t> & leb);

#endif
