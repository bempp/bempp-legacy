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
// and has been relicensed with permission by the HyENA authors. This does not affect the license of any
// other part of HyENA.


#ifndef galerkin_duffy_h
#define galerkin_duffy_h




// own includes
#include "singularitytraits.H"
#include "enumerators.H"
#include "mat.hpp"




/**
 * @ingroup galerkin
 */
template<ELEMENT_SHAPE SHAPE, SING_INT SINGULARITY>
class GalerkinDuffyExpression
{

  enum { num_regions = SingularityTraits<SHAPE, SINGULARITY>::num_gal_regions };

  typedef Point2  point_type[num_regions];
  typedef double weight_type[num_regions];


public:

  //! ctor
  template<typename QUADRATURE>
  GalerkinDuffyExpression(const QUADRATURE& qrule)
    :	num_quad_points(qrule.getNumPoints()),
      reference_points_x(NULL),
      reference_points_y(NULL),
      weights(NULL)
  {
    reference_points_x = new  point_type* [num_quad_points];
    reference_points_y = new  point_type* [num_quad_points];
    weights            = new weight_type* [num_quad_points];

    for (unsigned int q=0; q<num_quad_points; ++q) {
      reference_points_x[q] = new  point_type [num_quad_points];
      reference_points_y[q] = new  point_type [num_quad_points];
      weights[           q] = new weight_type [num_quad_points];
    }

    for (unsigned int i=0; i<num_quad_points; ++i)
      for (unsigned int j=0; j<num_quad_points; ++j) {
        this -> transform( reference_points_x[i][j],
                           reference_points_y[i][j],
                           weights[           i][j],
                           qrule.getPoint(i),
                           qrule.getPoint(j) );
      }
  }


  //! dtor
  ~GalerkinDuffyExpression()
  {
    for (unsigned int q=0; q<num_quad_points; ++q) {
      delete[] reference_points_x[q];
      delete[] reference_points_y[q];
      delete[] weights[           q];
    }

    delete[] reference_points_x;
    delete[] reference_points_y;
    delete[] weights;
  }






  const Point2& getPointX(const unsigned int i,
                          const unsigned int j,
                          const unsigned int k) const
  {
    assert(i<num_quad_points);
    assert(j<num_quad_points);
    assert(k<num_regions);
    return reference_points_x[i][j][k];
  }


  const Point2& getPointY(const unsigned int i,
                          const unsigned int j,
                          const unsigned int k) const
  {
    assert(i<num_quad_points);
    assert(j<num_quad_points);
    assert(k<num_regions);
    return reference_points_y[i][j][k];
  }


  const double getWeight(const unsigned int i,
                         const unsigned int j,
                         const unsigned int k) const
  {
    assert(i<num_quad_points);
    assert(j<num_quad_points);
    assert(k<num_regions);
    return weights[i][j][k];
  }


  const unsigned int getNumRegions() const
  {
    return num_regions;
  }




private:
  const unsigned int num_quad_points;

  point_type** reference_points_x;
  point_type** reference_points_y;

  weight_type** weights;


  /**
   * Expression of gauss points by using sauter expression rues
   */
  void transform(point_type& rpts_x,
                 point_type& rpts_y,
                 weight_type& wgts,
                 const Point2& x,
                 const Point2& y);
};















//=========================================================================
//
// TRIANGLE - TRIANGLE
//
//=========================================================================


// RE
template<>
void GalerkinDuffyExpression<TRIANGLE, REGULAR>::
transform (point_type& reference_points_x,
           point_type& reference_points_y,
           weight_type& weights,
           const Point2& x,
           const Point2& y)
{
  // some references to make things more readable
  const double & xsi  = x[ 0 ];
  const double & eta1 = x[ 1 ];
  const double & eta2 = y[ 0 ];
  const double & eta3 = y[ 1 ];

  // REGION 1
  reference_points_x[ 0 ][ 0 ] = xsi;
  reference_points_x[ 0 ][ 1 ] = xsi * eta1;
  reference_points_y[ 0 ][ 0 ] = eta2;
  reference_points_y[ 0 ][ 1 ] = eta2 * eta3;

  //------------------
  // pre-factors for each region
  weights[ 0 ] = xsi * eta2;
}


// CO
template<>
void GalerkinDuffyExpression<TRIANGLE, COINCIDENT>::
transform (point_type& reference_points_x,
           point_type& reference_points_y,
           weight_type& weights,
           const Point2& x,
           const Point2& y)
{
  // some references to make things more readable
  const double & xsi  = x[ 0 ];
  const double & eta1 = x[ 1 ];
  const double & eta2 = y[ 0 ];
  const double & eta3 = y[ 1 ];

  // some auxillary values
  const double eta123 = eta1 * eta2 * eta3;
  const double eta12  = eta1 * eta2;

  // REGION 1
  reference_points_x[ 0 ][ 0 ] = xsi;
  reference_points_x[ 0 ][ 1 ] = xsi * ( 1. - eta1 + eta12 );
  reference_points_y[ 0 ][ 0 ] = xsi * ( 1. - eta123 );
  reference_points_y[ 0 ][ 1 ] = xsi * ( 1. - eta1   );

  // REGION 2
  reference_points_x[ 1 ][ 0 ] = xsi * ( 1. - eta123 );
  reference_points_x[ 1 ][ 1 ] = xsi * ( 1. - eta1   );
  reference_points_y[ 1 ][ 0 ] = xsi;
  reference_points_y[ 1 ][ 1 ] = xsi * ( 1. - eta1 + eta12 );

  // REGION 3
  reference_points_x[ 2 ][ 0 ] = xsi;
  reference_points_x[ 2 ][ 1 ] = xsi * ( eta1  - eta12 + eta123 );
  reference_points_y[ 2 ][ 0 ] = xsi * (    1. - eta12 );
  reference_points_y[ 2 ][ 1 ] = xsi * ( eta1  - eta12 );

  // REGION 4
  reference_points_x[ 3 ][ 0 ] = xsi * (    1. - eta12 );
  reference_points_x[ 3 ][ 1 ] = xsi * ( eta1  - eta12 );
  reference_points_y[ 3 ][ 0 ] = xsi;
  reference_points_y[ 3 ][ 1 ] = xsi * ( eta1  - eta12 + eta123 );

  // REGION 5
  reference_points_x[ 4 ][ 0 ] = xsi * (    1. - eta123 );
  reference_points_x[ 4 ][ 1 ] = xsi * ( eta1  - eta123 );
  reference_points_y[ 4 ][ 0 ] = xsi;
  reference_points_y[ 4 ][ 1 ] = xsi * ( eta1  - eta12  );

  // REGION 6
  reference_points_x[ 5 ][ 0 ] = xsi;
  reference_points_x[ 5 ][ 1 ] = xsi * ( eta1  - eta12  );
  reference_points_y[ 5 ][ 0 ] = xsi * (    1. - eta123 );
  reference_points_y[ 5 ][ 1 ] = xsi * ( eta1  - eta123 );

  //------------------
  // pre-factors for each region
  weights[ 0 ] = xsi * xsi * xsi * eta1 * eta1 * eta2;
  weights[ 1 ] = weights[ 0 ];
  weights[ 2 ] = weights[ 0 ];
  weights[ 3 ] = weights[ 0 ];
  weights[ 4 ] = weights[ 0 ];
  weights[ 5 ] = weights[ 0 ];
}


// EA
template<>
void GalerkinDuffyExpression<TRIANGLE, EDGE_ADJACENT>::
transform (point_type& reference_points_x,
           point_type& reference_points_y,
           weight_type& weights,
           const Point2& x,
           const Point2& y)
{
  // some references to make things more readable
  const double & xsi  = x[ 0 ];
  const double & eta1 = x[ 1 ];
  const double & eta2 = y[ 0 ];
  const double & eta3 = y[ 1 ];
  //-------------------------------------------------------------
  // some auxilary values
  const double eta12  = eta1 * eta2;
  const double eta123 = eta1 * eta2 * eta3;

  // REGION 1
  reference_points_x[ 0 ][ 0 ] = xsi;
  reference_points_x[ 0 ][ 1 ] = xsi * eta1 * eta3;
  reference_points_y[ 0 ][ 0 ] = xsi * (    1. - eta12 );
  reference_points_y[ 0 ][ 1 ] = xsi * ( eta1  - eta12 );

  // REGION 2
  reference_points_x[ 1 ][ 0 ] = xsi;
  reference_points_x[ 1 ][ 1 ] = xsi * eta1;
  reference_points_y[ 1 ][ 0 ] = xsi * (     1. - eta123 );
  reference_points_y[ 1 ][ 1 ] = xsi * ( eta12  - eta123 );

  // REGION 3
  reference_points_x[ 2 ][ 0 ] = xsi * (    1. - eta12 );
  reference_points_x[ 2 ][ 1 ] = xsi * ( eta1  - eta12 );
  reference_points_y[ 2 ][ 0 ] = xsi;
  reference_points_y[ 2 ][ 1 ] = xsi * eta123;

  // REGION 4
  reference_points_x[ 3 ][ 0 ] = xsi * (     1. - eta123 );
  reference_points_x[ 3 ][ 1 ] = xsi * ( eta12  - eta123 );
  reference_points_y[ 3 ][ 0 ] = xsi;
  reference_points_y[ 3 ][ 1 ] = xsi * eta1;

  // REGION 5
  reference_points_x[ 4 ][ 0 ] = xsi * (    1. - eta123 );
  reference_points_x[ 4 ][ 1 ] = xsi * ( eta1  - eta123 );
  reference_points_y[ 4 ][ 0 ] = xsi;
  reference_points_y[ 4 ][ 1 ] = xsi * eta12;

  //------------------
  // pre-factors for each region
  weights[ 0 ] = xsi * xsi * xsi * eta1 * eta1;
  weights[ 1 ] = weights[ 0 ] * eta2;
  weights[ 2 ] = weights[ 1 ];
  weights[ 3 ] = weights[ 1 ];
  weights[ 4 ] = weights[ 1 ];
}




// VA
template<>
void GalerkinDuffyExpression<TRIANGLE, VRTX_ADJACENT>::
transform (point_type& reference_points_x,
           point_type& reference_points_y,
           weight_type& weights,
           const Point2& x,
           const Point2& y)
{
  // some references to make things more readable
  const double & xsi  = x[ 0 ];
  const double & eta1 = x[ 1 ];
  const double & eta2 = y[ 0 ];
  const double & eta3 = y[ 1 ];

  // REGION 1
  reference_points_x[ 0 ][ 0 ] = xsi;
  reference_points_x[ 0 ][ 1 ] = xsi * eta1;
  reference_points_y[ 0 ][ 0 ] = xsi * eta2;
  reference_points_y[ 0 ][ 1 ] = xsi * eta2 * eta3;

  // REGION 2
  reference_points_x[ 1 ][ 0 ] = xsi * eta2;
  reference_points_x[ 1 ][ 1 ] = xsi * eta2 * eta3;
  reference_points_y[ 1 ][ 0 ] = xsi;
  reference_points_y[ 1 ][ 1 ] = xsi * eta1;

  //------------------
  // pre-factors for each region
  weights[ 0 ] = xsi * xsi * xsi * eta2;
  weights[ 1 ] = weights[ 0 ];
}






//=========================================================================
//
// QUADRANGLE - QUADRANGLE
//
//=========================================================================


// RE
template<>
void GalerkinDuffyExpression<QUADRANGLE, REGULAR>::
transform (point_type& reference_points_x,
           point_type& reference_points_y,
           weight_type& weights,
           const Point2& x,
           const Point2& y)
{
  // some references to make things more readable
  const double & xsi  = x[ 0 ];
  const double & eta1 = x[ 1 ];
  const double & eta2 = y[ 0 ];
  const double & eta3 = y[ 1 ];

  // REGION 1
  reference_points_x[ 0 ][ 0 ] = xsi;
  reference_points_x[ 0 ][ 1 ] = eta1;
  reference_points_y[ 0 ][ 0 ] = eta2;
  reference_points_y[ 0 ][ 1 ] = eta3;

  //------------------
  // pre-weights for each region
  weights[ 0 ] = 1.;
}


// CO
template<>
void GalerkinDuffyExpression<QUADRANGLE, COINCIDENT>::
transform (point_type& reference_points_x,
           point_type& reference_points_y,
           weight_type& weights,
           const Point2& x,
           const Point2& y)
{
  // some references to make things more readable
  const double & xsi  = x[ 0 ];
  const double & eta1 = x[ 1 ];
  const double & eta2 = y[ 0 ];
  const double & eta3 = y[ 1 ];

  // some auxiliary data
  const double aux1 = ( 1. - xsi ) * eta3;
  const double aux2 = ( 1. - xsi * eta1 ) * eta2;
  const double aux3 = xsi + aux1;
  const double aux4 = xsi * eta1 + aux2;

  // REGION 1
  reference_points_x[ 0 ][ 0 ] = aux1;
  reference_points_x[ 0 ][ 1 ] = aux2;
  reference_points_y[ 0 ][ 0 ] = aux3;
  reference_points_y[ 0 ][ 1 ] = aux4;

  // REGION 2
  reference_points_x[ 1 ][ 0 ] = aux2;
  reference_points_x[ 1 ][ 1 ] = aux1;
  reference_points_y[ 1 ][ 0 ] = aux4;
  reference_points_y[ 1 ][ 1 ] = aux3;

  // REGION 3
  reference_points_x[ 2 ][ 0 ] = aux1;
  reference_points_x[ 2 ][ 1 ] = aux4;
  reference_points_y[ 2 ][ 0 ] = aux3;
  reference_points_y[ 2 ][ 1 ] = aux2;

  // REGION 4
  reference_points_x[ 3 ][ 0 ] = aux2;
  reference_points_x[ 3 ][ 1 ] = aux3;
  reference_points_y[ 3 ][ 0 ] = aux4;
  reference_points_y[ 3 ][ 1 ] = aux1;

  // REGION 5
  reference_points_x[ 4 ][ 0 ] = aux3;
  reference_points_x[ 4 ][ 1 ] = aux2;
  reference_points_y[ 4 ][ 0 ] = aux1;
  reference_points_y[ 4 ][ 1 ] = aux4;

  // REGION 6
  reference_points_x[ 5 ][ 0 ] = aux4;
  reference_points_x[ 5 ][ 1 ] = aux1;
  reference_points_y[ 5 ][ 0 ] = aux2;
  reference_points_y[ 5 ][ 1 ] = aux3;

  // REGION 7
  reference_points_x[ 6 ][ 0 ] = aux3;
  reference_points_x[ 6 ][ 1 ] = aux4;
  reference_points_y[ 6 ][ 0 ] = aux1;
  reference_points_y[ 6 ][ 1 ] = aux2;

  // REGION 8
  reference_points_x[ 7 ][ 0 ] = aux4;
  reference_points_x[ 7 ][ 1 ] = aux3;
  reference_points_y[ 7 ][ 0 ] = aux2;
  reference_points_y[ 7 ][ 1 ] = aux1;

  //------------------
  // pre-weights for each region
  weights[ 0 ] = xsi * ( 1. - xsi ) * ( 1. - xsi * eta1 );
  weights[ 1 ] = weights[ 0 ];
  weights[ 2 ] = weights[ 0 ];
  weights[ 3 ] = weights[ 0 ];
  weights[ 4 ] = weights[ 0 ];
  weights[ 5 ] = weights[ 0 ];
  weights[ 6 ] = weights[ 0 ];
  weights[ 7 ] = weights[ 0 ];
}


// EA
template<>
void GalerkinDuffyExpression<QUADRANGLE, EDGE_ADJACENT>::
transform (point_type& reference_points_x,
           point_type& reference_points_y,
           weight_type& weights,
           const Point2& x,
           const Point2& y)
{
  // some references to make things more readable
  const double & xsi  = x[ 0 ];
  const double & eta1 = x[ 1 ];
  const double & eta2 = y[ 0 ];
  const double & eta3 = y[ 1 ];

  // REGION 1
  reference_points_x[ 0 ][ 0 ] = ( 1. - xsi ) * eta3 + xsi;
  reference_points_x[ 0 ][ 1 ] = xsi * eta2;
  reference_points_y[ 0 ][ 0 ] = ( 1. - xsi ) * eta3;
  reference_points_y[ 0 ][ 1 ] = xsi * eta1;

  // REGION 2
  reference_points_x[ 1 ][ 0 ] = ( 1. - xsi ) * eta3;
  reference_points_x[ 1 ][ 1 ] = xsi * eta2;
  reference_points_y[ 1 ][ 0 ] = xsi + ( 1. - xsi ) * eta3;
  reference_points_y[ 1 ][ 1 ] = xsi * eta1;

  // REGION 3
  reference_points_x[ 2 ][ 0 ] = ( 1. - xsi * eta1 ) * eta3 + xsi * eta1;
  reference_points_x[ 2 ][ 1 ] = xsi * eta2;
  reference_points_y[ 2 ][ 0 ] = ( 1. - xsi * eta1 ) * eta3;
  reference_points_y[ 2 ][ 1 ] = xsi;

  // REGION 4
  reference_points_x[ 3 ][ 0 ] = ( 1. - xsi * eta1 ) * eta3 + xsi * eta1;
  reference_points_x[ 3 ][ 1 ] = xsi;
  reference_points_y[ 3 ][ 0 ] = ( 1. - xsi * eta1 ) * eta3;
  reference_points_y[ 3 ][ 1 ] = xsi * eta2;

  // REGION 5
  reference_points_x[ 4 ][ 0 ] = ( 1. - xsi * eta1 ) * eta3;
  reference_points_x[ 4 ][ 1 ] = xsi * eta2;
  reference_points_y[ 4 ][ 0 ] = ( 1. - xsi * eta1 ) * eta3 + xsi * eta1;
  reference_points_y[ 4 ][ 1 ] = xsi;

  // REGION 6
  reference_points_x[ 5 ][ 0 ] = ( 1. - xsi * eta1 ) * eta3;
  reference_points_x[ 5 ][ 1 ] = xsi;
  reference_points_y[ 5 ][ 0 ] = ( 1. - xsi * eta1 ) * eta3 + xsi * eta1;
  reference_points_y[ 5 ][ 1 ] = xsi * eta2;

  //------------------
  // pre-weights for each region
  weights[ 0 ] = xsi * xsi * ( 1. - xsi );
  weights[ 1 ] = weights[ 0 ];
  weights[ 2 ] = xsi * xsi * ( 1. - xsi * eta1 );
  weights[ 3 ] = weights[ 2 ];
  weights[ 4 ] = weights[ 2 ];
  weights[ 5 ] = weights[ 2 ];
}



// VA
template<>
void GalerkinDuffyExpression<QUADRANGLE, VRTX_ADJACENT>::
transform (point_type& reference_points_x,
           point_type& reference_points_y,
           weight_type& weights,
           const Point2& x,
           const Point2& y)
{

  // some references to make things more readable
  const double & xsi  = x[ 0 ];
  const double & eta1 = x[ 1 ];
  const double & eta2 = y[ 0 ];
  const double & eta3 = y[ 1 ];

  const double aux1 = xsi * eta1;
  const double aux2 = xsi * eta2;
  const double aux3 = xsi * eta3;

  // REGION 1
  reference_points_x[ 0 ][ 0 ] = xsi;
  reference_points_x[ 0 ][ 1 ] = aux1;
  reference_points_y[ 0 ][ 0 ] = aux2;
  reference_points_y[ 0 ][ 1 ] = aux3;

  // REGION 2
  reference_points_x[ 1 ][ 0 ] = aux1;
  reference_points_x[ 1 ][ 1 ] = xsi;
  reference_points_y[ 1 ][ 0 ] = aux2;
  reference_points_y[ 1 ][ 1 ] = aux3;

  // REGION 3
  reference_points_x[ 2 ][ 0 ] = aux1;
  reference_points_x[ 2 ][ 1 ] = aux2;
  reference_points_y[ 2 ][ 0 ] = xsi;
  reference_points_y[ 2 ][ 1 ] = aux3;

  // REGION 4
  reference_points_x[ 3 ][ 0 ] = aux1;
  reference_points_x[ 3 ][ 1 ] = aux2;
  reference_points_y[ 3 ][ 0 ] = aux3;
  reference_points_y[ 3 ][ 1 ] = xsi;

  //------------------
  // pre-weights for each region
  weights[ 0 ] = xsi * xsi * xsi;
  weights[ 1 ] = weights[ 0 ];
  weights[ 2 ] = weights[ 0 ];
  weights[ 3 ] = weights[ 0 ];

}




//////////////////////////////////////////////////////////////////////
#endif //galerkin_duffy_h
