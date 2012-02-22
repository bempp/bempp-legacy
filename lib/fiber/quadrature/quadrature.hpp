// Copyright (C) 2009-2010 Matthias Messner, Michael Messner, Franz
// Rammerstorfer, Peter Urthaler
// 
// This file is part of HyENA - a C++ boundary element methods library.
// 
// HyENA is free software: you can redistribute it and/or modify it under the
// terms of the GNU Lesser Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
// 
// HyENA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser Public License for more
// details.
// 
// You should have received a copy of the GNU Lesser Public License along with
// HyENA. If not, see <http://www.gnu.org/licenses/>.

/**
 * @file   quadrature.hpp
 * @author Michael
 * @date   created:     26.08.09
 *         last change: 11.12.09
 */
#ifndef quadrature_hpp
#define quadrature_hpp





// own includeds  
#include "gausstensor.hpp"
#include "gausstria.hpp"
#include "gauss1D.hpp"


//#include "macros.H"
//#include "enumerators.H"
#include "traits.H"
#include "mat.hpp"






/**
 * @ingroup quadrature
 *
 * A QuadraturePoint provides its coordinates on the reference element and
 * the related quadrature weight.
 * @tparam NUM_COORDS
 */
template<int NUM_COORDS>
struct QuadraturePoint
{
	double coords[NUM_COORDS];
	double quadrature_weight;
};



/**
 * @ingroup quadrature
 *
 * The QuadratureRule object is a non specialzed template class for all
 * different @p QUADRATURE_RULES and @p ELEMENT_SHAPES (e.g. GAUSS for
 * TRIANGLES). 
 * @tparam ELEMENT_SHAPE 
 * @tparam QUADRATURE_RULE
 */
template<ELEMENT_SHAPE SHAPE, QUADRATURE_RULE RULE>
class QuadratureRule
{
	/**
	 * Dimension of the problem is defined by ELEMENT_SHAPE, which is known at
	 * compile time.
	 */
  enum{ shape_dim = ShapeTraits<SHAPE>::shape_dim };



	/**
	 * Definition of local coordinate system trough the dimension of the
	 * problem at compile time.
	 */
	typedef typename PointTraits<shape_dim>::point_type      local_point_type;









public:
	/**
	 * This constructor (quadrature order known) is to be specialzed for
	 * different template parameters.  
	 * @param[in] order quadrature order
	 */
	QuadratureRule(unsigned int order);



	/**
	 * Own destructor to free quad_point_array
	 */
	~QuadratureRule( )
	{ }





	/**
	 * The function getPoint returns the local coordinates for a chosen 
	 * QuadraturePoint.
	 * @param[in] node 
	 * @return local_point_type local coordinates of QuadraturePoint
	 */
	const local_point_type getPoint(unsigned int node) const;

 

	/**
	 * The function getWeight returns the quadrature weight for a chosen 
	 * QuadraturePoint.
	 * @param[in] node 
	 * @return double quadrature weight
	 */
	const double getWeight(unsigned int node) const;



	/**
	 * The function getNumPoints returns the  number of quadrature points per 
	 * Element.
	 * @return unsigned int number of QuadraturePoint 's
	 */
	const unsigned int getNumPoints() const
  { return num_of_points; }






	const local_point_type getPoint(unsigned int order,unsigned int node) const;
	const double getWeight(unsigned int order, unsigned int node) const;
	const unsigned int getNumPoints(unsigned int order) const;





private:	
	/**
	 * Standard constructor not to be used
	 */
	QuadratureRule( );


	/**
	 * The copy constructor is intentionally not defined in order to prevent
	 * QuadratureRule from being copied.
	 */
 	QuadratureRule(const QuadratureRule&);



	/**
	 * The assignment operator is intentionally not defined because must not be
	 * used.
	 */
	const QuadratureRule& operator= (const QuadratureRule&);




	unsigned int order;


	unsigned int num_of_points;


	unsigned int position;

};
























//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////                                         ////////////////
/////////////      S P E C I A L I Z A T  I O N       ////////////////
/////////////                                         ////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Specializations for LINES.


/////////////////////////////////////////////////////////////////////
// Specializations for TRIANGLES.

template<> inline
QuadratureRule<TRIANGLE, GAUSS>::QuadratureRule(unsigned int _order)
	: order(_order), 
		num_of_points(gTriaPointsPerOrder[order-1]),
		position(gTriaAddress[gTriaPointsPerOrder[order-1]-1])
{ }


template<> inline
const Point2 QuadratureRule<TRIANGLE, GAUSS>::getPoint(unsigned int node) const
{
	return Point2(gpTria[(position+ node)*3+1]+
								gpTria[(position+ node)*3+2],
								gpTria[(position+ node)*3+2]);
}


template<> inline
const double QuadratureRule<TRIANGLE, GAUSS>::getWeight(unsigned int node) const
{
	return 0.5*gwTria[position+ node];
}



template<> inline
const Point2 QuadratureRule<TRIANGLE, GAUSS>::getPoint(unsigned int order,
																											 unsigned int node) const
{
	unsigned int pos( (gTriaAddress[gTriaPointsPerOrder[order-1]-1] + node)*3 );
	return Point2(gpTria[pos+1]+
								gpTria[pos+2],
								gpTria[pos+2]);
}

template<> inline
const double QuadratureRule<TRIANGLE, GAUSS>::getWeight(unsigned int order,
																												unsigned int node) const
{
	return 0.5*gwTria[gTriaAddress[gTriaPointsPerOrder[order-1]-1]+node];
}


template<> inline
const unsigned int QuadratureRule<TRIANGLE, GAUSS>::
getNumPoints(unsigned int order) const
{
	return gTriaPointsPerOrder[order-1];
}

//////////////////////////////////////////////////////////////////////
// Specializations for QUADRANGLES.

template<> inline
QuadratureRule<QUADRANGLE, GAUSS>::QuadratureRule(unsigned int _order)
	: order(_order), 
		num_of_points(order*order),
		position(gTensorAddress[order-1])
{ }




template<> inline
const Point2 QuadratureRule<QUADRANGLE, GAUSS>::
		getPoint(unsigned int node) const
{
	return Point2(gpTensor[(position+ node)*2  ],
								gpTensor[(position+ node)*2+1]);
}


template<> inline
const double QuadratureRule<QUADRANGLE, GAUSS>::
getWeight(unsigned int node) const
{
	return gwTensor[position+node];
}

template<> inline
const Point2 QuadratureRule<QUADRANGLE, GAUSS>::
getPoint(unsigned int order, unsigned int node) const
{
	unsigned int pos( (gTensorAddress[order-1]+ node)*2 );
	
	return Point2(gpTensor[pos  ],
								gpTensor[pos+1]);
}


template<> inline
const double QuadratureRule<QUADRANGLE, GAUSS>::
getWeight(unsigned int order, unsigned int node) const
{
	return gwTensor[gTensorAddress[order-1]+ node];
}

template<> inline
const unsigned int QuadratureRule<QUADRANGLE, GAUSS>::
getNumPoints(unsigned int order) const
{
	return order*order;
}



/// **
//  * Specialized constructor for HYPERCUBE.
//  * @tparam HYPERCUBE
//  * @tparam GAUSS
//  * @param[in] order quadrature order
//  */
// template<> inline
// void QuadratureRule<HYPERCUBE,GAUSS>::setUpQuadrature(unsigned int order)
// {
//  QuadratureRule<QUADRANGLE, GAUSS> quad;

//  num_of_points = quad.getNumPoints()*quad.getNumPoints();
//  quad_point_array = new QuadraturePoint<4>[num_of_points];

// 	for(unsigned int i=0; i<quad.getNumPoints(); ++i)
//   {
// 	  for(unsigned int j=0; j<quad.getNumPoints(); ++j)
//     {
//       (quad_point_array[j+quad.getNumPoints()*i]).coords[0]=quad.getPoint(i)[0];
//       (quad_point_array[j+quad.getNumPoints()*i]).coords[1]=quad.getPoint(i)[1];
//       (quad_point_array[j+quad.getNumPoints()*i]).coords[2]=quad.getPoint(j)[0];
//       (quad_point_array[j+quad.getNumPoints()*i]).coords[3]=quad.getPoint(j)[1];
		
// 		  quad_point_array[j+quad.getNumPoints()*i].quadrature_weight= 
// 				quad.getWeight(i)*quad.getWeight(j);
//     }
// 	}		
// }



#endif
