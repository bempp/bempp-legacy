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
 * @file   scalartraits.H
 * @ingroup traits
 *
 * every type exept specialisations isn't a scalar-type!
 * The @p ScalarTraits return true if the type @p T is a scalar.
 * Feel free to add your type!  
 *
 * @tparam T type
 *
 * @author Rf
 * @date   created:     02.10.09
 *         last change: 09.12.09
 */
#ifndef scalartraits_h
#define scalartraits_h



// system includes
#include <complex>
#include <cmath>



// base struct
template<typename T> struct ScalarTraits
{	enum{is_scalar = false}; };

template<> 
struct ScalarTraits<int>
{	enum{is_scalar = true}; };

template<> 
struct ScalarTraits<float>
{	enum{is_scalar = true}; };

template<> 
struct ScalarTraits<double>
{	enum{is_scalar = true}; };

template<> 
struct ScalarTraits<std::complex<float> >
{	enum{is_scalar = true}; };

template<> 
struct ScalarTraits<std::complex<double> >
{	enum{is_scalar = true}; };



#endif // include guard

