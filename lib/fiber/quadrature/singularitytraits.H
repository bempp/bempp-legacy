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
 * @file   singularitytraits.H
 * @ingroup traits
 *
 * The @p SingularityTraits - class defines compile time constants for
 * specific singularity types
 *
 * @tparam ELEMENT_SHAPE element shape
 * @tparam SING_INT      singularity type
 *
 * @author Matthias, Michael
 * @date   created:     02.02.10
 *         last change: 02.02.10
 */
#ifndef singularitytraits_h
#define singularitytraits_h



// own include
#include "enumerators.H"
#include "shapetraits.H"



template <ELEMENT_SHAPE SHAPE, SING_INT SINGULARITY>
struct SingularityTraits;


template<ELEMENT_SHAPE SHAPE> struct SingularityTraits<SHAPE, REGULAR>
{
	enum {num_gal_regions = ShapeTraits<SHAPE>::gal_re_regions};
};


template<ELEMENT_SHAPE SHAPE> struct SingularityTraits<SHAPE, VRTX_ADJACENT>
{
	enum {num_col_regions = ShapeTraits<SHAPE>::col_va_regions,
				num_gal_regions = ShapeTraits<SHAPE>::gal_va_regions};
};

template<ELEMENT_SHAPE SHAPE> struct SingularityTraits<SHAPE, EDGE_ADJACENT>
{
	enum {num_col_regions = ShapeTraits<SHAPE>::col_ea_regions,
				num_gal_regions = ShapeTraits<SHAPE>::gal_ea_regions};
};

template<ELEMENT_SHAPE SHAPE> struct SingularityTraits<SHAPE, COINCIDENT>
{
	enum {num_col_regions = ShapeTraits<SHAPE>::col_co_regions,
				num_gal_regions = ShapeTraits<SHAPE>::gal_co_regions};
};  




#endif
