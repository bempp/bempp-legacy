// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "max_distance.hpp"

#include "grid.hpp"

#include <cmath>

namespace Bempp
{

double maxDistance(const Grid& grid1, const Grid& grid2)
{
    const int dimWorld = grid1.dimWorld();
    if (grid2.dimWorld() != dimWorld)
        throw std::invalid_argument("maxDistance(): both grids must be "
                                    "embedded in spaces of the same dimension");

    arma::Col<double> lowerBound1, upperBound1;
    grid1.getBoundingBox(lowerBound1, upperBound1);
    arma::Col<double> lowerBound2, upperBound2;
    grid2.getBoundingBox(lowerBound2, upperBound2);

    if (dimWorld == 3) {
        const int pointCount = 16; // 2 * 2**dimWorld
        arma::Mat<double> points(dimWorld, pointCount);
        points.col(0)  = lowerBound1;
        points.col(1)  = lowerBound1; points(0, 1) = upperBound1(0);
        points.col(2)  = lowerBound1; points(1, 2) = upperBound1(1);
        points.col(3)  = lowerBound1; points(2, 3) = upperBound1(2);
        points.col(4)  = upperBound1; points(0, 4) = lowerBound1(0);
        points.col(5)  = upperBound1; points(1, 5) = lowerBound1(1);
        points.col(6)  = upperBound1; points(2, 6) = lowerBound1(2);
        points.col(7)  = upperBound1;
        points.col(8)  = lowerBound2;
        points.col(9)  = lowerBound2; points(0,  9) = upperBound2(0);
        points.col(10) = lowerBound2; points(1, 10) = upperBound2(1);
        points.col(11) = lowerBound2; points(2, 11) = upperBound2(2);
        points.col(12) = upperBound2; points(0, 12) = lowerBound2(0);
        points.col(13) = upperBound2; points(1, 13) = lowerBound2(1);
        points.col(14) = upperBound2; points(2, 14) = lowerBound2(2);
        points.col(15) = upperBound2;

        double result = 0.;
        for (int i = 0; i < pointCount; ++i)
            for (int j = i + 1; j < pointCount; ++j) {
                arma::Col<double> diff = points.col(i) - points.col(j);
                double dist = sqrt(diff(0) * diff(0) + diff(1) * diff(1) +
                                   diff(2) * diff(2));
                result = std::max(dist, result);
            }

        return result;
    }
    else
        throw std::runtime_error("maxDistance(): currently implemented only "
                                 "for grids embedded in 3D spaces");
}

} // namespace Bempp
