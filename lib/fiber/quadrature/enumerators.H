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

//file: enumerators.H
//--------------------------------------------------------------------------
#ifndef enumerators_h
#define enumerators_h
//--------------------------------------------------------------------------
/** \file
 *  Contains Enumerators for the BEM code
 */
//--------------------------------------------------------------------------

/** Enumerator for element shapes */
enum ELEMENT_SHAPE{
    /* 1 */ NO_SHAPE,       //!< element shape not (yet) set
    /* 1 */ LINE,           //!< line element
	/* 2 */ TRIANGLE,       //!< triangle element
    /* 3 */ QUADRANGLE      //!< quadrangle element
};


/** Enumerator describing the Type of Singular Integration */
enum SING_INT {
	/* 0 */ NO_SING_INT,    //!< singularity type is not set yet
	/* 1 */ REGULAR,        //!< Regular integration
	/* 2 */ COINCIDENT,     //!< coincident
	/* 3 */ EDGE_ADJACENT,  //!< edge adjacent
	/* 4 */ VRTX_ADJACENT  //!< vertex adjacent
};

/** Enumerator for quadrature rule */
enum QUADRATURE_RULE {
	/* 0 */ NO_RULE,          //!< Quadrature Rule not (yet) set
    /* 1 */ GAUSS             //!<
};

#endif // enumerators_h
