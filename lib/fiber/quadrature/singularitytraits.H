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
