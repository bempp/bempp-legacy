// Copyright (C) 2009-2010 Matthias Messner, Michael Messner, Franz
// Rammerstorfer, Peter Urthaler
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

// Note: This file is originally part of the HyENA project
// (http://portal.tugraz.at/portal/page/portal/Files/i2610/files/Forschung/Software/HyENA/html/index.html)
// and has been relicensed with permission by the HyENA authors. This does not
// affect the license of any
// other part of HyENA.

/**
 * @file   shapetraits.H
 * @ingroup traits
 *
 * The @p ShapeTraits - class defines compile time constants for specific
 * element shapes
 *
 * @tparam ELEMENT_SHAPE element shape
 *
 * @author Matthias, Michael, Rf
 * @date   created:     01.08.09
 *         last change: 09.12.09
 */
#ifndef shapetraits_h
#define shapetraits_h

// system include
#include <complex>

// own include
#include "enumerators.H"

// base struct
template <ELEMENT_SHAPE SHAPE> struct ShapeTraits;

// LINE
template <> struct ShapeTraits<LINE> {
  enum {
    dim = 2,            // dimension of domain
    shape_dim = 1,      // dimension of reference boundary element
    num_coords = 1,     // number coordinates
    num_topo_nodes = 2, // number of topology nodes
    identity_scale = 2, // scale identity
    num_edges = 2
  };

  // CollocationPoints on DISCONTINUOUS boundaries need to be indented
  static constexpr double indent = 0.25;
};

// TRIANGLE
template <> struct ShapeTraits<TRIANGLE> {
  enum {
    dim = 3,            // dimension of domain
    shape_dim = 2,      // dimension of reference boundary element
    num_coords = 3,     // barycentric coordinates
    num_topo_nodes = 3, // number of topology nodes
    identity_scale = 2, // scale identity
    gal_re_regions = 1, // integration: sauter regions: regular
    gal_va_regions = 2, // integration: sauter regions: vertex adjacent
    gal_ea_regions = 5, // integration: sauter regions: edge adjacent
    gal_co_regions = 6, // integration: sauter regions: coincident
    col_va_regions = 1, // integration: collo regions: vertex adjacent
    col_ea_regions = 2, // integration: collo regions: edge adjacent
    col_co_regions = 3, // integration: collo regions: coincident
    num_edges = 3
  };

  // CollocationPoints on DISCONTINUOUS boundaries need to be indented
  static constexpr double indent = 0.16666666666;
};

// QUADRANGLE
template <> struct ShapeTraits<QUADRANGLE> {
  enum {
    dim = 3,            // dimension of domain
    shape_dim = 2,      // dimension of reference boundary element
    num_coords = 2,     // number coordinates
    num_topo_nodes = 4, // number of topology nodes
    identity_scale = 1, // scale identity
    gal_re_regions = 1, // integration: sauter regions: regular
    gal_va_regions = 4, // integration: sauter regions: vertex adjacent
    gal_ea_regions = 6, // integration: sauter regions: edge adjacent
    gal_co_regions = 8, // integration: sauter regions: coincident
    col_va_regions = 2, // integration: collo regions: vertex adjacent
    col_ea_regions = 4, // integration: collo regions: edge adjacent
    col_co_regions = 4, // integration: collo regions: coincident
    num_edges = 4
  };

  // CollocationPoints on DISCONTINUOUS boundaries need to be indented
  static constexpr double indent = 0.125;
};

//// HYPERCUBE
// template<> struct ShapeTraits<HYPERCUBE>
//{
//  enum
//		{dim            = 3,  // dimension of domain
//		 shape_dim      = 4,  // dimension of reference boundary element
//		 num_coords     = 4,  // number coordinates
//		 num_topo_nodes = 8,  // number of topology nodes
//		 identity_scale = 1}; // scale identity
//};

#endif // include guard
