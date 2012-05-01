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
 * @file   mat.hpp
 * @ingroup linalg
 * @author Rf
 * @date   created:     08.09.09
 *         last change: 11.11.09
 */
#ifndef mat_hpp
#define mat_hpp



// system includes
#include <iostream>
#include <cmath>
#include <complex>



// own includes
#include "macros.H"
#include "traits.H"



/**
 * Provide a dense matrix class with compile-time fixed length.
 * Since the length is known this @p Mat is more efficient than
 * matrices with runtime-dependent length.
 * storage type: row major
 *
 * This class is also used for column-vectors of any value type and for
 * geometric points, which define a location in dim-dimensional real space.
 * For this purpose global typedefs are implemented: @p Point2, @p Point3,
 * Vec4d and so on.
 * Point_ is the preferred object to be passed to functions, don't use C-arrays.
 *
 * Note: @p T aka. value_type can be of Mat-type too (or any other type),
 * so you can create nested matrices/vectors. ATTENTION: Most functions are
 * implemented for matices with scalar entries only, therefore nested
 * matrices are just like a container.
 *
 * @tparam T value_type
 * @tparam NUM_ROWS number of rows
 * @tparam NUM_COLS number of colums
 */
template<typename T, int NUM_ROWS, int NUM_COLS>
class Mat
{
public:


  //! @name compile time constants
  //@{
  typedef T value_type;                    //!< type of entries
  enum {
    num_rows   = NUM_ROWS,                 //!< number of rows
    num_cols   = NUM_COLS,                 //!< number of columns
    size       = NUM_ROWS*NUM_COLS,        //!< total number of entries
    is_vector  = ( NUM_COLS == 1 ),        //!< column-vector
    is_square  = ( NUM_ROWS == NUM_COLS ), //!< square matrix
  };
	//@}



  //////////////////////////////////////////////////////////////////////
  //   CONSTRUCTORS / DESTRUCTOR   /////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  /**
   * standard constructor. no initialisation of values.
   */
  Mat( ){/*empty*/}



  /**
   * ctor with value
   * @param[in] _val set all values of mat to @p _val
   *
   * e.g. use this constructor to init Mat with zero:
   * Mat<double,3,2> dummy( NumberTraits<value_type>::zero() );
   */
  explicit Mat( const value_type _val)
  {
    for(unsigned int i=0; i<size; ++i)
      vals[i] = _val;
  }



  /**
   * ctor with C-array
   * @param[in] _vals C-array of values
   */
  explicit Mat( const value_type _vals[size] )
  {
    for(unsigned int i=0; i<size; ++i)
      vals[i] = _vals[i];
  }



  /**
   * copy ctor
   * @param[in] m Mat to be copied
   */
  Mat(const Mat& m)
  {
    for(unsigned int i=0; i<size; ++i)
      vals[i] = m.vals[i];
  }



  /**
   * ctor for 2D-point
   * @param[in] x,y coordinates
   */
  explicit Mat(const double x,
      const double y);



  /**
   * ctor for 3D-point
   * @param[in] x,y,z coordinates
   */
  explicit Mat(const double x,
      const double y,
      const double z);



  /**
   * standard dtor
   */
  ~Mat( ){/*empty*/}



  //////////////////////////////////////////////////////////////////////
  //   OVERLOADING OPERATORS   /////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  /**
   * const acess to components. (read only)
   * @param[in] i,j row index, column index
   * @return const i-j-th value
   */
  const value_type& operator () (const unsigned int i,
                                 const unsigned int j) const
  {
    HYENA_ASSERT( i < num_rows );
    HYENA_ASSERT( j < num_cols );
    return vals[ j*num_rows+i ];
  }



  /**
   * acess to components. (read/write)
   * @param[in] i,j row index, column index
   * @return i-j-th value
   */
  value_type& operator () (const unsigned int i,
                           const unsigned int j)
  {
    HYENA_ASSERT( i < num_rows );
    HYENA_ASSERT( j < num_cols );
    return vals[ j*num_rows+i ];
  }



  /**
   * const acess to components intended for vectors. (read only)
   * @param[in] id storage index (row major)
   * @return const id-th value
   */
  const value_type& operator [] (const unsigned int id) const
  {
    HYENA_ASSERT( id < size );
    return vals[id];
  }



  /**
   * acess to components intended for vectors. (read/write)
   * @param[in] id storage index (row major)
   * @return const id-th value
   */
  value_type& operator [] (const unsigned int id)
  {
    HYENA_ASSERT( id < size );
    return vals[id];
  }



  /**
   * assignment operator
   * @param[in] m mat to be assigned
   * @return modified mat
   */
  Mat& operator = (const Mat& m)
  {
    for(unsigned int i=0; i<size; ++i)
      vals[i] = m.vals[i];
    return *this;
  }



  /**
   * assignment operator
   * @tparam TYPE scalar-type of value @p s. e.g. int, double,...
   * @param[in] s scalar-value assigned to each entry of Mat.
   * @return modified mat
   */
  template <typename TYPE>
  Mat& operator = (const TYPE& s)
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    for(unsigned int i=0; i<size; ++i)
      vals[i] = s;
    return *this;
  }


  /**
   * compare two matrices
   * @param[in] m matrix to be compared
   * @return equality (true/false)
   */
  bool operator == (const Mat& m)	const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    unsigned int counter =0;
    for(unsigned int i=0; i<size; ++i){
      double diff = NumberTraits<value_type>::abs( vals[i]-m.vals[i] );
      if ( diff <= NumberTraits<value_type>::tolerance() )
        counter++;
    }
    if(counter == size)
      return true;
    else
      return false;
  }



  /**
   * compare two matrices
   * @param[in] m matrix to be compared
   * @return inequality (true/false)
   */
  bool operator != (const Mat& m)	const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    unsigned int counter =0;
    for(unsigned int i=0; i<size; ++i){
      double diff = NumberTraits<value_type>::abs( vals[i]-m.vals[i] );
      if ( diff <= NumberTraits<value_type>::tolerance() )
        counter++;
    }
    if(counter != size)
      return true;
    else
      return false;
  }



  /**
   * add one mat to the other
   * @param[in] m mat to be added
   * @return modified mat
   */
  Mat& operator += (const Mat& m)
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    for(unsigned int i=0; i<size; ++i)
      vals[i] += m.vals[i];
    return *this;
  }



  /**
   * subtract one mat from the other
   * @param[in] m mat to be subtracted
   * @return modified mat
   */
  Mat& operator -= (const Mat& m)
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    for(unsigned int i=0; i<size; ++i)
      vals[i] -= m.vals[i];
    return *this;
  }



  /**
   * scalar multiplication
   * @tparam TYPE scalar-type of value @p s. e.g. int, double,...
   * @param[in] s scalar-value assigned to each entry of Mat.
   * @return modified mat
   * NOTE: only post-multiplications of type "Mat * scalar" are possible
   */
  template<typename TYPE>
  Mat& operator *= (const TYPE s)
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    for(unsigned int i=0; i<size; ++i)
      vals[i] *= s;
    return *this;
  }



  /**
   * scalar division
   * @tparam TYPE scalar-type of value @p s. e.g. int, double,...
   * @param[in] s scalar-value assigned to each entry of Mat.
   * @return modified mat
   * NOTE: only post-divisions of type "Mat * scalar" are possible
   */
  template<typename TYPE>
  Mat& operator /= (const TYPE s)
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    for(unsigned int i=0; i<size; ++i)
      vals[i] /= s;
    return *this;
  }



  /**
   * + operator: add two mats
   * @param[in] m mat to be added to
   * @return new mat
   */
  Mat operator + (const Mat& m) const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    value_type res[size];
    for(unsigned int i=0; i<size; ++i)
      res[i] = vals[i] + m.vals[i];
    // create object with array constructor
    return Mat(res);
  }



  /**
   * - operator: subtract two mats
   * @param[in] m mat to be subtracted of
   * @return new mat
   */
  Mat operator - (const Mat& m) const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    value_type res[size];
    for(unsigned int i=0; i<size; ++i)
      res[i] = vals[i] - m.vals[i];
    // create object with array constructor
    return Mat(res);
  }



  /**
   * Unary minus operator. Negate all
   * entries of a mat.
   * @return new mat
   */
  Mat operator - () const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    value_type res[size];
    for(unsigned int i=0; i<size; ++i)
      res[i] = -vals[i];
    // create object with array constructor
    return Mat(res);
  }



  /**
   * scalar multiplication
   * @tparam TYPE scalar-type of value @p s. e.g. int, double,...
   * @param[in] s scalar-value assigned to each entry of Mat.
   * @return new mat
   * NOTE: only post-multiplications of type "Mat * scalar" are possible
   */
  template<typename TYPE>
  Mat operator * (const TYPE s)	const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    value_type res[size];
    for(unsigned int i=0; i<size; ++i)
      res[i] = vals[i] * s;
    // create object with array constructor
    return Mat(res);
  }



  /**
   * scalar division
   * @tparam TYPE scalar-type of value @p s. e.g. int, double,...
   * @param[in] s scalar-value assigned to each entry of Mat.
   * @return new mat
   * NOTE: only post-divisions of type "Mat / scalar" are possible
   */
  template<typename TYPE>
  Mat operator / (const TYPE s)	const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    value_type res[size];
    for(unsigned int i=0; i<size; ++i)
      res[i] = vals[i] / s;
    // create object with array constructor
    return Mat(res);
  }



  /**
   * universal matrix-matrix multiplication
   * \f$ \mathbf{C}=\mathbf{A} \cdot \mathbf{B}\f$ or
   * \f$ C_{ik}= \sum\limits_j A_{ij}B_{jk}\f$
   * aplicable for matrix-vector multiplication
   * @tparam ROWS rows of matrix m
   * @tparam COLS columns of matrix m
   * @param[in] m matrix to be multiplied
   * @return new Mat
   */
  template <int ROWS,int COLS>
  Mat<value_type,num_rows,COLS> operator * (const Mat<value_type,
                                            ROWS,COLS>& m) const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    HYENA_STATIC_ASSERT(ROWS == num_cols,BAD_INDEX);
    Mat<value_type,num_rows,COLS> res;
    for(unsigned i=0; i<num_rows; ++i)
      for(unsigned int j=0; j<num_cols; ++j)
        for(unsigned int k=0; k<COLS; ++k)
          res(i,k) += (*this)(i,j)*m(j,k);
    return res;
  }



  //////////////////////////////////////////////////////////////////////
  //   FUNCTIONS   /////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  // universal functions
  //--------------------------------------------------------------------


  //! set all values to zero.
  void zeros()
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    for(unsigned int i=0; i<size; ++i)
      vals[i] = NumberTraits<value_type>::zero();
  }



  //! set all values to one.
  void ones()
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    for(unsigned int i=0; i<size; ++i)
      vals[i] = NumberTraits<value_type>::one();
  }



  //! @return number of rows
  const int getNumRows() const {return num_rows;}



  //! @return number of columns
  const int getNumCols() const {return num_cols;}



  //! @return total number of entries (size=num_rows*num_columns)
  const int getNumEntries() const {return size;}




  //! write out matrix in octave format
  std::ostream& writeOctaveFormat(std::ostream& out = std::cout) const
  {
    if(is_vector) {
      for (unsigned int i=0; i<size; ++i)
        out << vals[i] << std::endl;
      return out;

    }
    else {
      for (unsigned int row=0; row<num_rows; ++row) {
        for (unsigned int col=0; col<num_cols; ++col)
          out << (*this)(row, col) << " ";
        out << std::endl;
      }
    }
    return out;
  }



  // functions for rectangular matices
  //--------------------------------------------------------------------

  /**
   * transpose
   * use this function only for non-quadratic matrices!
   * @return new mat
   */
  Mat<value_type,num_cols,num_rows> trans() const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    Mat<value_type,num_cols,num_rows> res;
    for(unsigned int i=0; i<num_rows; ++i)
      for(unsigned int j=0; j<num_cols; ++j)
        res(j,i) = (*this)(i,j);
    return res;
  }



  // functions for square matices
  //--------------------------------------------------------------------

  //! build identity matrix.
  void identity()
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    HYENA_STATIC_ASSERT_SQUARE_MATRIX(Mat);
    this->zeros();
    for(unsigned int i=0; i<num_rows; ++i)
      vals[i*(1+num_rows)] = NumberTraits<value_type>::one();
  }



  /**
   * transpose for square matrices.
   * store result A <- A^T
   */
  void qtrans()
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    HYENA_STATIC_ASSERT_SQUARE_MATRIX(Mat);
    Mat<value_type,num_rows,num_cols> res;
    for(unsigned int i=0; i<num_rows; ++i)
      for(unsigned int j=0; j<num_cols; ++j)
        res(j,i) = (*this)(i,j);
    (*this)=res;
  }



  /**
   * matrix-matrix multiplication for square matrices.
   * \f$ \mathbf{C}=\mathbf{A} \cdot \mathbf{B}\f$ or
   * \f$ C_{ik}= \sum\limits_j A_{ij}B_{jk}\f$
   * store result B <- A.B
   * @param[in,out] m matrix to be multiplied
   */
  void qmultM(Mat& m)
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    HYENA_STATIC_ASSERT_SQUARE_MATRIX(Mat);
    Mat<value_type,num_rows,num_cols> res;
    for(unsigned i=0; i<num_rows; ++i)
      for(unsigned int j=0; j<num_cols; ++j)
        for(unsigned int k=0; k<num_rows; ++k)
          res(i,k) += (*this)(i,j)*m(j,k);
    m=res;
  }



  /**
   * matrix-vector multiplication for square matrices.
   * \f$ \mathbf{C}=\mathbf{A} \cdot \mathbf{B}\f$ or
   * \f$ C_{ik}= \sum\limits_j A_{ij}B_{jk}\f$
   * store result b <- A.b
   * @param[in,out] v vector to be multiplied
   */
  void qmultV(Mat<value_type,num_rows,1>& v)
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    HYENA_STATIC_ASSERT_SQUARE_MATRIX(Mat);
    Mat<value_type,num_rows,1> res;
    for(unsigned i=0; i<num_rows; ++i)
      for(unsigned int j=0; j<num_cols; ++j)
        res[i] += (*this)(i,j)*v[j];
    v=res;
  }



  // functions for vectors
  //--------------------------------------------------------------------

  /**
   * dot (inner) product of two vectors.
   * \f$ a = <\mathbf{u},\mathbf{v}> = {\sum\limits_i u_i v_i} \f$
   * @param[in] v vector to be multiplied
   * @return value of norm
   */
  value_type dotProduct(const Mat& v) const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    HYENA_STATIC_ASSERT_VECTOR(Mat);
    value_type dot = NumberTraits<value_type>::zero();
    for(unsigned int i=0; i<size; ++i)
      dot += vals[i] * v.vals[i];
    return dot;
  }



  /**
   * dyadic (outer) product of two vectors.
   * \f$ \mathbf{A}=\mathbf{u} \otimes \mathbf{v}\f$ or
   * \f$ A_{ij} =  u_i v_j \f$.
   * @param[in] v vector to be multiplied
   * @return new Mat
   */
  Mat<value_type,num_rows,num_rows> dyadicProduct(const Mat& v) const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    HYENA_STATIC_ASSERT_VECTOR(Mat);
    Mat<value_type,num_rows,num_rows> res;
    for(unsigned i=0; i<num_rows; ++i)
      for(unsigned int j=0; j<num_rows; ++j)
        res(i,j) = vals[i] * v.vals[j];
    return res;
  }


  /**
   * \f$ L_2\f$ -norm of a vector.
   * \f$ |\mathbf{v}|_2 = \sqrt{\sum\limits_i v_i^2} \f$
   * @return value of norm
   */
  double norm() const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    HYENA_STATIC_ASSERT_VECTOR(Mat);
    double norm2 = 0.;
    for(unsigned int i=0; i<size; ++i)
      norm2 += NumberTraits<value_type>::abs( vals[i] )
        * NumberTraits<value_type>::abs( vals[i] );
    return sqrt( norm2 );
  }



  /**
   * Returns the Euclidian distance of
   * <tt>this</tt> to <tt>p</tt>, i.e. the <tt>L_2</tt> norm
   * of the difference vector between the two geometric points.
   * NOTE: this function works only for real valued vectors.
   * @param[in] p second point
   * @return value of distance
   */
  double distance(const Mat& p) const
  {
    HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
    HYENA_STATIC_ASSERT_VECTOR(Mat);
    HYENA_STATIC_ASSERT_REAL_NUMBER(value_type);
    double sum = 0.;
    for (unsigned int i=0; i<size; ++i) {
      const double diff = vals[i]-p.vals[i];
      sum += diff*diff;
    }
    return sqrt(sum);
  }



  /**
   * cross product between two vectors.
   * specialisation only for Mat<double,3,1> aka. Point3
   * @param[in] v second vector
   * @return new vector
   */
  Mat crossProduct(const Mat& v) const;



private:

  value_type vals[size]; //!< values, stored in a simple C-array

};



//////////////////////////////////////////////////////////////////////
//  SPECIALISATION SECTION   /////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------
// constructor for 2d-point
template<> inline
Mat<double,2,1>::Mat(const double x,
                     const double y)
{
  vals[0] = x;
  vals[1] = y;
}

//--------------------------------------------------------------------
// constructor for 3d-point
template<> inline
Mat<double,3,1>::Mat(const double x,
                     const double y,
                     const double z)
{
  vals[0] = x;
  vals[1] = y;
  vals[2] = z;
}



//--------------------------------------------------------------------
// cross product for T == double and size ==3 aka. Point3 only!
template<> inline
Mat<double,3,1> Mat<double,3,1>::crossProduct(const Mat& v) const
{
  value_type _vals[] = {0.,0.,0.};
  // compute the product
  _vals[0] = vals[1]*v.vals[2] - vals[2]*v.vals[1];
  _vals[1] = vals[2]*v.vals[0] - vals[0]*v.vals[2];
  _vals[2] = vals[0]*v.vals[1] - vals[1]*v.vals[0];
  return Mat<double,3,1>(_vals);
}



//////////////////////////////////////////////////////////////////////
//   GLOBAL TYPEDEFS   ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


/**
 * @name Global typedefs for geometric points,etc.
 *
 * Note that the typedefs have following structure:
 * appended number: specifies size e.g. @p Mat3 for 3x3 matrix
 * appended character: specifies value_type e.g. cd for complex<double>
 */
//@{
// points
typedef Mat<double,1,1>               Point1;
typedef Mat<double,2,1>               Point2;
typedef Mat<double,3,1>               Point3;
typedef Mat<double,4,1>               Point4;
// vectors int
typedef Mat<int,1,1>                  Vec1i;
typedef Mat<int,2,1>                  Vec2i;
typedef Mat<int,3,1>                  Vec3i;
typedef Mat<int,4,1>                  Vec4i;
// vectors double
typedef Mat<double,1,1>               Vec1d;
typedef Mat<double,2,1>               Vec2d;
typedef Mat<double,3,1>               Vec3d;
typedef Mat<double,4,1>               Vec4d;
typedef Mat<double,6,1>               Vec6d;
typedef Mat<double,9,1>               Vec9d;
// vectors complex<double>
typedef Mat<std::complex<double>,2,1> Vec2cd;
typedef Mat<std::complex<double>,3,1> Vec3cd;
typedef Mat<std::complex<double>,4,1> Vec4cd;
//square matices int
typedef Mat<int,1,1>                  Mat1i;
typedef Mat<int,2,2>                  Mat2i;
typedef Mat<int,3,3>                  Mat3i;
typedef Mat<int,4,4>                  Mat4i;
//square matices double
typedef Mat<double,1,1>               Mat1d;
typedef Mat<double,2,2>               Mat2d;
typedef Mat<double,3,3>               Mat3d;
typedef Mat<double,4,4>               Mat4d;
//square matices complex<double>
typedef Mat<std::complex<double>,1,1> Mat1cd;
typedef Mat<std::complex<double>,2,2> Mat2cd;
typedef Mat<std::complex<double>,3,3> Mat3cd;
typedef Mat<std::complex<double>,4,4> Mat4cd;
//@}



//////////////////////////////////////////////////////////////////////
//   POINT TRAITS   //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/**
 * base struct
 * @tparam DIM dimension of point
 */
template<int DIM> struct PointTraits;



/**
 * @name PointTraits
 *
 * The @p PointTraits return the point type dependend on the space dimension.
 */
//@{
template<>
struct PointTraits<1>
{	typedef Mat<double,1,1> point_type; };

template<>
struct PointTraits<2>
{	typedef Mat<double,2,1> point_type; };

template<>
struct PointTraits<3>
{	typedef Mat<double,3,1> point_type; };

template<>
struct PointTraits<4>
{	typedef Mat<double,4,1> point_type; };

template<>
struct PointTraits<5>
{	typedef Mat<double,5,1> point_type; };

template<>
struct PointTraits<6>
{	typedef Mat<double,6,1> point_type; };
//@}


//////////////////////////////////////////////////////////////////////
//   OUTPUT   ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/**
 * overloading output operator.
 * @tparam T value_type
 * @tparam NUM_ROWS number of rows
 * @tparam NUM_COLS number of colums
 * @param[out] out output-stream
 * @param[in] m mat to write
 */
template<typename T, int NUM_ROWS, int NUM_COLS>
std::ostream& operator << (std::ostream& out,
                           const Mat<T,NUM_ROWS,NUM_COLS>& m)
{
  if( m.is_vector )	{
    out << '[' << NUM_ROWS*NUM_COLS << "](";
    if (NUM_ROWS*NUM_COLS > 0)
      out << m[0];
    for(unsigned int i=1; i<NUM_ROWS*NUM_COLS; ++i)
      out << ','<< m[i];
    out << ')';
  }
  else {
    out << '[' << NUM_ROWS << ',' << NUM_COLS << "](";
    if (NUM_ROWS > 0) {
      out << '(' ;
      if (NUM_COLS > 0)
        out << m(0, 0);
      for (unsigned int j = 1; j < NUM_COLS; ++ j)
        out << ',' << m(0, j);
      out << ')';
    }
    for (unsigned int i = 1; i < NUM_ROWS; ++ i) {
      out << ",(" ;
      if (NUM_COLS > 0)
        out << m(i, 0);
      for (unsigned int j = 1; j < NUM_COLS; ++ j)
        out << ',' << m(i, j);
      out << ')';
    }
    out << ')';
  }
  return out;
}

template <typename TYPE, typename T,int NUM_ROWS,int NUM_COLS>
Mat<T,NUM_ROWS,NUM_COLS> operator* (TYPE s,const Mat<T,NUM_ROWS,NUM_COLS>& m)
{
  typedef typename Mat<T,NUM_ROWS,NUM_COLS>::value_type value_type;
  HYENA_STATIC_ASSERT_SCALAR_TYPE(value_type);
  unsigned size = Mat<T,NUM_ROWS,NUM_COLS>::size;
  value_type res[size];
  for(unsigned int i=0; i<size; ++i)
    res[i] = m[i] * s;
  // create object with array constructor
  return Mat<T,NUM_ROWS,NUM_COLS>(res);
}



#endif //include guard
