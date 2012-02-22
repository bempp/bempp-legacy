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
