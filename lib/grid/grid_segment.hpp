// Copyright (C) 2011-2013 by the BEM++ Authors
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

#ifndef bempp_grid_segment_hpp
#define bempp_grid_segment_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"

#include <boost/array.hpp>
#include <set>
#include <vector>

namespace Bempp
{

class Grid;

class GridSegment
{
public:
    static GridSegment wholeGrid(const Grid& grid);
    static GridSegment openDomain(const Grid& grid, int domain);
    static GridSegment closedDomain(const Grid& grid, int domain);

    GridSegment(const Grid& grid,
                const std::set<int>& excludedEntitiesCodim0,
                const std::set<int>& excludedEntitiesCodim1,
                const std::set<int>& excludedEntitiesCodim2,
                const std::set<int>& excludedEntitiesCodim3);

    GridSegment(int entityCountCodim0,
                int entityCountCodim1,
                int entityCountCodim2,
                int entityCountCodim3,
                const std::set<int>& excludedEntitiesCodim0,
                const std::set<int>& excludedEntitiesCodim1,
                const std::set<int>& excludedEntitiesCodim2,
                const std::set<int>& excludedEntitiesCodim3);

    bool contains(int codim, int index) const;

    void markExcludedEntities(int codim, std::vector<int>& marks,
                              int mark = -1) const;

    GridSegment complement() const;
    GridSegment union_(const GridSegment& other) const;
    GridSegment difference(const GridSegment& other) const;
    GridSegment intersection(const GridSegment& other) const;

//    Segment operator+(const Segment& other) const;
//    Segment operator-(const Segment& other) const;
//    Segment operator*(const Segment& other) const;

private:
    boost::array<int, 4> m_entityCounts;
    boost::array<std::set<int>, 4> m_excludedEntities;
};

GridSegment gridSegmentWithPositiveX(const Grid& grid);

} // namespace Bempp

#endif
