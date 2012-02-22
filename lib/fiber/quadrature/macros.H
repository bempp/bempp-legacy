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
 * @file   macros.H
 *         Provide some useful macros that behave like functions 
 *         in the code and define exeptions.
 * @ingroup exceptions
 *
 * @author Rf
 * @date   created:     10.09.09
 *         last change: 10.09.09
 */
#ifndef macros_h
#define macros_h



// system includes
#include <stdexcept>

//own includes
#include "numbertraits.H"
#include "scalartraits.H"

//////////////////////////////////////////////////////////////////////
//   OpenMP defines  /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

#ifdef STATIC_COLLO_OMP
#define ASSEMBLE_COLLO_OMP
#endif

#ifdef STATIC_GALERKIN_OMP
#define ASSEMBLE_GALERKIN_OMP
#endif

#ifdef SAUTER_OMP
#define FREQUENCY_LOOP_OMP
#define FUNDSOL_OMP
#endif

#ifdef LUBICH_COLLO_OMP
#define ASSEMBLE_COLLO_OMP
#define FUNDSOL_OMP
#endif

#ifdef LUBICH_GALERKIN_OMP
#define ASSEMBLE_GALERKIN_OMP
#define FUNDSOL_OMP
#endif



//////////////////////////////////////////////////////////////////////
//   EXCEPTIONS   ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//! @exception  HyenaLogicError logic error exeption (standard case)
class HyenaLogicError : public std::logic_error
{
public:
  HyenaLogicError(const std::string& msg) 
    : std::logic_error(msg) {}
};



//! @exception HyenaErrorMsg logic error with message 
class HyenaErrorMsg : public std::logic_error
{
public:
  HyenaErrorMsg(const std::string& msg) 
    : std::logic_error(msg) {}
};



//! @exception HyenaNotImplemented not(jet) implemented
class HyenaNotImplemented : public std::logic_error
{
public:
  HyenaNotImplemented(const std::string& msg) 
    : std::logic_error(msg) {}
};



//! @exception HyenaFileError file doesn't exist or isn't readable
class HyenaFileError : public std::runtime_error
{
public:
  HyenaFileError(const std::string& msg) 
    : std::runtime_error(msg) {}
};



//! throw exeption
#define HYENA_THROW(e) do { throw e; } while (0)





//////////////////////////////////////////////////////////////////////
//   RUNTIME ASSERTION   /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/**
 * The HYENA_ASSERT() macro acts like C's assert(), but throws a 
 * HYENA_ERROR() instead of just exiting.
 * Note that this marco acts only in DEBUG mode.
 */
#ifdef DEBUG
#define HYENA_ASSERT(asserted)                                          \
  do {if (!(asserted)){                                                 \
      BOOST_LOG_TRIVIAL(error) << "Assertion '" #asserted "' failed." ; \
      HYENA_ERROR(); } } while(0)
#else
#define HYENA_ASSERT(asserted) //do nothing in RELEASE mode.
#endif

/**
 * If you want to make sure you are accessing a section of code just
 * stick a HYENA_HERE(); in it.
 */
#define HYENA_HERE()                                            \
  do { BOOST_LOG_TRIVIAL(debug) << __FILE__ << ", line "        \
                                << __LINE__ << ", compiled "    \
                                << __DATE__ << " at "           \
                                << __TIME__ ; } while (0)



/**
 * The HYENA_ERROR() macro throws a  HyenaLogicError exception
 */
#define HYENA_ERROR()                                   \
  do { HYENA_HERE();                                    \
    const std::string msg = "internal logic";           \
    BOOST_LOG_TRIVIAL(error) << msg;                    \
    HYENA_THROW(HyenaLogicError(msg)); } while(0)



/**
 * The HYENA_ERROR_MSG(const std::string& msg) macro throws a HyenaErrorMsg
 * exception.
 */
#define HYENA_ERROR_MSG(msg)                    \
  do { HYENA_HERE();                            \
    BOOST_LOG_TRIVIAL(error) << msg;            \
    HYENA_THROW(HyenaErrorMsg(msg)); } while(0)


								 
/**
 * The HYENA_NOT_IMPLEMENTED() macro throws a HyenaNotImplemented exception.
 */
#define HYENA_NOT_IMPLEMENTED()                         \
  do { HYENA_HERE();                                    \
    const std::string msg = "not implemented";          \
    BOOST_LOG_TRIVIAL(error) << msg;                    \
    HYENA_THROW(HyenaNotImplemented(msg)); } while(0)



/**
 * The HYENA_FILE_ERROR(const std::string& filename) macro throws a
 * HyenaFileError exception.
 */
#define HYENA_FILE_ERROR(filename)                              \
  do { HYENA_HERE();                                            \
    const std::string msg = "cannot access file: "+ filename;   \
    BOOST_LOG_TRIVIAL(error) << msg;                            \
    HYENA_THROW(HyenaFileError(msg)); } while(0)



 


//////////////////////////////////////////////////////////////////////
//   COMPILE TIME ASSERTION   ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/**
 * This generates a compile time error message if the condition is not true.
 * In other words it is the compile time equivalent of the assert macro.
 * Note that if the condition is true, then the macro will generate 
 * neither code nor data. When used in a template, the static assertion 
 * will be evaluated at the time the template is instantiated, 
 * therefore it is particularly useful for validating template parameters. 
 * Note that this marco acts in DEBUG and RELEASE mode.
 *
 * @tparam CONDITION must be a compile time boolean expression
 */
template<bool CONDITION>
struct hyena_static_assert {};

//! specialisation for true
template<>
struct hyena_static_assert<true>
{
  enum {
    ONLY_VECTORS,
    ONLY_SQUARE_MATRICES,
    ONLY_SCALAR_TYPES,
    ONLY_REAL_NUMBERS,
    ONLY_3D,
    BAD_INDEX,
    WRONG_MATERIAL_PARAMETER,
    DIFFERENT_ELEMENTS,
    YOU_MADE_A_PROGRAMMING_MISTAKE,
    YOU_POOR_FOOL,
    NOT_YET_IMPLEMENTED
  };
};

/**
 * static assertion.
 * @param[in] CONDITION to be asserted
 * @param[in] MSG an enum listed in struct hyena_static_assert<true>
 *
 * Example:
 * \code
 * HYENA_STATIC_ASSERT(dim==2, YOU_MADE_A_PROGRAMMING_MISTAKE);
 * \endcode
 *
 * feel free to add your personal HYENA_STATIC_ASSERT_... and your
 * individual error message! 
 */
#define HYENA_STATIC_ASSERT(CONDITION,MSG)                      \
  if (hyena_static_assert<CONDITION ? true : false>::MSG) {}

//! static assertion failing if the @p TYPE is not a vector type
#define HYENA_STATIC_ASSERT_VECTOR(TYPE)        \
  HYENA_STATIC_ASSERT(TYPE::is_vector,          \
                      ONLY_VECTORS)

//! static assertion failing if the @p TYPE is not a square matrix type
#define HYENA_STATIC_ASSERT_SQUARE_MATRIX(TYPE) \
  HYENA_STATIC_ASSERT(TYPE::is_square,          \
                      ONLY_SQUARE_MATRICES)

//! static assertion failing if the @p TYPE is not a real number
#define HYENA_STATIC_ASSERT_REAL_NUMBER(TYPE)                   \
  HYENA_STATIC_ASSERT(NumberTraits<TYPE>::is_complex==false,    \
                      ONLY_REAL_NUMBERS)

//! static assertion failing if the @p TYPE is not a scalar type
#define HYENA_STATIC_ASSERT_SCALAR_TYPE(TYPE)           \
  HYENA_STATIC_ASSERT(ScalarTraits<TYPE>::is_scalar,    \
                      ONLY_SCALAR_TYPES)



#endif // include guard
