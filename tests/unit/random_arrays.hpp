// Copyright (C) 2011 by the BEM++ Authors
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

#ifndef bempp_generate_random_matrix_hpp
#define bempp_generate_random_matrix_hpp

#include "common/armadillo_fwd.hpp"
#include <boost/type_traits/is_complex.hpp>
#include <boost/utility/enable_if.hpp>
#include "common/eigen_support.hpp"

using namespace Bempp;

// Real ValueType
template <typename ValueType>
typename boost::disable_if<
boost::is_complex<ValueType>,
Vector<ValueType>
>::type
generateRandomVector(int rowCount)
{
    return Vector<ValueType>::Random(rowCount);
}

// Complex ValueType
template <typename ValueType>
typename boost::enable_if<
boost::is_complex<ValueType>,
Vector<ValueType>
>::type
generateRandomVector(int rowCount)
{
    return Vector<ValueType>::Random(rowCount) +
            ValueType(0, 1) *
            Vector<ValueType>::Random(rowCount);
}

// Real ValueType
template <typename ValueType>
typename boost::disable_if<
boost::is_complex<ValueType>,
Matrix<ValueType>
>::type
generateRandomMatrix(int rowCount, int colCount)
{
    return Matrix<ValueType>::Random(rowCount, colCount);
}

// Complex ValueType
template <typename ValueType>
typename boost::enable_if<
boost::is_complex<ValueType>,
Matrix<ValueType>
>::type
generateRandomMatrix(int rowCount, int colCount)
{
    return Matrix<ValueType>::Random(rowCount, colCount) +
            ValueType(0, 1) *
            Matrix<ValueType>::Random(rowCount, colCount);
}

#endif
