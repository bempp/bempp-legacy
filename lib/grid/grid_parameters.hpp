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

#ifndef bempp_grid_parameters_hpp
#define bempp_grid_parameters_hpp

namespace Bempp {

/** \ingroup grid
    \brief %Grid parameters.

  This structure is used to specify parameters of grid constructed by GridFactory.
  */
struct GridParameters {
    /** \brief %Grid topology */
    enum Topology {
        /** \brief one-dimensional grid embedded in a two-dimensional space */
        LINEAR,
        /** \brief two-dimensional grid composed of triangular elements,
            embedded in a three-dimensional space */
        TRIANGULAR,
        /** \brief two-dimensional grid composed of quadrilateral elements,
            embedded in a three-dimensional space */
        QUADRILATERAL,
        /** \brief two-dimensional grid composed (potentially) both of
            triangular and quadrilateral elements, embedded in a
            three-dimensional space */
        HYBRID_2D,
        /** \brief three-dimensional grid composed of tetrahedral elements,
            embedded in a three-dimensional space*/
        TETRAHEDRAL
    } topology;
};

} // namespace Bempp

#endif
