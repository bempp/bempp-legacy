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
template<ELEMENT_SHAPE SHAPE> struct ShapeTraits;



// LINE
template<> struct ShapeTraits<LINE>
{
	enum
		{dim            = 2,  // dimension of domain                     
		 shape_dim      = 1,  // dimension of reference boundary element 
		 num_coords     = 1,  // number coordinates                      
		 num_topo_nodes = 2,  // number of topology nodes                
		 identity_scale = 2, // scale identity
	   num_edges = 2};

	// CollocationPoints on DISCONTINUOUS boundaries need to be indented
	static const double indent = 0.25;
};



// TRIANGLE
template<> struct ShapeTraits<TRIANGLE>
{
	enum
		{dim            = 3,  // dimension of domain
		 shape_dim      = 2,  // dimension of reference boundary element 
		 num_coords     = 3,  // barycentric coordinates
		 num_topo_nodes = 3,  // number of topology nodes
		 identity_scale = 2,  // scale identity
		 gal_re_regions = 1,  // integration: sauter regions: regular         
		 gal_va_regions = 2,  // integration: sauter regions: vertex adjacent 
		 gal_ea_regions = 5,  // integration: sauter regions: edge adjacent   
		 gal_co_regions = 6,  // integration: sauter regions: coincident      
		 col_va_regions = 1,  // integration: collo regions: vertex adjacent  
		 col_ea_regions = 2,	// integration: collo regions: edge adjacent    
		 col_co_regions = 3, // integration: collo regions: coincident
	   num_edges = 3};
	
	// CollocationPoints on DISCONTINUOUS boundaries need to be indented
	static const double indent = 0.16666666666;
};



// QUADRANGLE
template<> struct ShapeTraits<QUADRANGLE>
{
	enum
		{dim            = 3,  // dimension of domain
		 shape_dim      = 2,  // dimension of reference boundary element 
		 num_coords     = 2,  // number coordinates
		 num_topo_nodes = 4,  // number of topology nodes
		 identity_scale = 1,  // scale identity 
		 gal_re_regions = 1,  // integration: sauter regions: regular          
		 gal_va_regions = 4,  // integration: sauter regions: vertex adjacent  
		 gal_ea_regions = 6,  // integration: sauter regions: edge adjacent    
		 gal_co_regions = 8,  // integration: sauter regions: coincident       
		 col_va_regions = 2,  // integration: collo regions: vertex adjacent  
		 col_ea_regions = 4,	// integration: collo regions: edge adjacent    
		 col_co_regions = 4, // integration: collo regions: coincident
	   num_edges = 4};
	
	// CollocationPoints on DISCONTINUOUS boundaries need to be indented
	static const double indent = 0.125;
};



//// HYPERCUBE
//template<> struct ShapeTraits<HYPERCUBE>
//{
//  enum
//		{dim            = 3,  // dimension of domain
//		 shape_dim      = 4,  // dimension of reference boundary element
//		 num_coords     = 4,  // number coordinates
//		 num_topo_nodes = 8,  // number of topology nodes
//		 identity_scale = 1}; // scale identity
//};



#endif // include guard
