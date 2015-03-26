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

#ifndef bempp_check_arrays_are_close_hpp
#define bempp_check_arrays_are_close_hpp

#include "fiber/scalar_traits.hpp"
#include "fiber/_2d_array.hpp"
#include "fiber/_3d_array.hpp"
#include "fiber/_4d_array.hpp"

#include "common/eigen_support.hpp"
#include <complex>
#include <iomanip>
#include <boost/test/unit_test.hpp>

using namespace Bempp;

template <typename ValueType>
boost::test_tools::predicate_result
check_arrays_are_close(const Matrix<ValueType>& left,
                       const Matrix<ValueType>& right,
                       typename Fiber::ScalarTraits<ValueType>::RealType tolerance)
{
    typedef typename Fiber::ScalarTraits<ValueType>::RealType RealType;
    const int digits10 = std::numeric_limits<RealType>::digits10;

    boost::test_tools::predicate_result result(true);
    if (left.rows() != right.rows() || left.cols() != right.cols()) {
        result = false;
        result.message() << "Size mismatch [("
                         << left.rows() << ", " << left.cols() << ") != ("
                         << right.rows() << ", " << right.cols() << ")]";
        return result;
    }
    for (size_t r = 0; r < left.rows(); ++r)
        for (size_t c = 0; c < left.cols(); ++c) {
            RealType diff = std::abs(left(r, c) - right(r, c));
            RealType avg = std::abs(left(r, c) + right(r, c)) / 2.;
            if (diff > tolerance * (1 + avg)) {
                result = false;
                result.message() << std::setprecision(digits10)
                                 << "\n  Mismatch at position ("
                                 << r << ", " << c << ") ["
                                 << left(r, c) << " != " << right(r, c) << "]";
            }
        }

    return result;
}

template <typename ValueType>
boost::test_tools::predicate_result
check_arrays_are_close(const std::vector<Matrix<ValueType>>& left,
                       const std::vector<Matrix<ValueType>>& right,
                       typename Fiber::ScalarTraits<ValueType>::RealType tolerance)
{
    typedef typename Fiber::ScalarTraits<ValueType>::RealType RealType;
    const int digits10 = std::numeric_limits<RealType>::digits10;

    boost::test_tools::predicate_result result(true);

    if (left.size()!=right.size()){
        result = false;
        result.message() << "Number of slices matrices mismatch: ["
                       <<  left.size() << " != " << right.size() << "]";
        return result;
    }

    for (size_t s = 0; s < left.size(); ++s){
        if (left[s].rows() != right[s].rows() ||
                left[s].cols() != right[s].cols()){
            result = false;
            result.message() << "Size mismatch in slice " << s << ". "
                           << "[(" << left[s].rows() << "," << left[s].cols() << ") != "
                           << "(" << right[s].rows() << "," << right[s].cols() << ")]";
        }
        for (size_t r = 0; r < left[s].rows(); ++r)
            for (size_t c = 0; c < left[s].cols(); ++c){
                RealType diff = std::abs(left[s](r,c) - right[s](r,c));
                RealType avg = std::abs(left[s](r,c) + right[s](r,c))/2.;
                if (diff > tolerance * (1+avg)) {
                    result = false;
                    result.message() << std::setprecision(digits10)
                                     << "\n  Mismatch in slice " << s << " at position ("
                                     << r << ", " << c << ") ["
                                     << left[s](r, c) << " != "
                                     << right[s](r, c) << "]";

                }
            }
    }

    return result;
}

template <typename ValueType>
boost::test_tools::predicate_result
check_arrays_are_close(const Fiber::_4dArray<ValueType>& left,
                       const Fiber::_4dArray<ValueType>& right,
                       typename Fiber::ScalarTraits<ValueType>::RealType tolerance)
{
    typedef typename Fiber::ScalarTraits<ValueType>::RealType RealType;
    const int digits10 = std::numeric_limits<RealType>::digits10;

    boost::test_tools::predicate_result result(true);
    if (left.extent(0) != right.extent(0) ||
            left.extent(1) != right.extent(1) ||
            left.extent(2) != right.extent(2) ||
            left.extent(3) != right.extent(3)) {
        result = false;
        result.message() << "Size mismatch [("
                         << left.extent(0) << ", " << left.extent(1) << ", "
                         << left.extent(2) << ", " << left.extent(3) << ") != ("
                         << right.extent(0) << ", " << right.extent(1) << ", "
                         << right.extent(2) << ", " << right.extent(3) << "]";
        return result;
    }
    for (size_t r = 0; r < left.extent(0); ++r)
        for (size_t c = 0; c < left.extent(1); ++c)
            for (size_t s = 0; s < left.extent(2); ++s)
                for (size_t u = 0; u < left.extent(3); ++u) {
                    RealType diff =
                            std::abs(left(r, c, s, u) - right(r, c, s, u));
                    RealType avg =
                            std::abs(left(r, c, s, u) + right(r, c, s, u)) / 2.;
                    if (diff > tolerance * (1 + avg)) {
                        result = false;
                        result.message() << std::setprecision(digits10)
                                         << "\n  Mismatch at position ("
                                         << r << ", " << c << ", "
                                         << s << ", " << u << ") ["
                                         << left(r, c, s, u) << " != "
                                         << right(r, c, s, u) << "]";
                    }
                }

    return result;
}

template <typename ValueType>
boost::test_tools::predicate_result
check_arrays_are_close(const Fiber::_3dArray<ValueType>& left,
                       const Fiber::_3dArray<ValueType>& right,
                       typename Fiber::ScalarTraits<ValueType>::RealType tolerance)
{
    typedef typename Fiber::ScalarTraits<ValueType>::RealType RealType;
    const int digits10 = std::numeric_limits<RealType>::digits10;

    boost::test_tools::predicate_result result(true);
    if (left.extent(0) != right.extent(0) ||
            left.extent(1) != right.extent(1) ||
            left.extent(2) != right.extent(2)) {
        result = false;
        result.message() << "Size mismatch [("
                         << left.extent(0) << ", " << left.extent(1) << ", "
                         << left.extent(2) << ") != ("
                         << right.extent(0) << ", " << right.extent(1) << ", "
                         << right.extent(2)  << "]";
        return result;
    }
    for (size_t r = 0; r < left.extent(0); ++r)
        for (size_t c = 0; c < left.extent(1); ++c)
            for (size_t s = 0; s < left.extent(2); ++s) {
                    RealType diff = std::abs(left(r, c, s) - right(r, c, s));
                    RealType avg = std::abs(left(r, c, s) + right(r, c, s)) / 2.;
                    if (diff > tolerance * (1 + avg)) {
                        result = false;
                        result.message() << std::setprecision(digits10)
                                         << "\n  Mismatch at position ("
                                         << r << ", " << c << ", "
                                         << s << ") ["
                                         << left(r, c, s) << " != "
                                         << right(r, c, s) << "]";
                    }
            }

    return result;
}

template <typename ValueType>
boost::test_tools::predicate_result
check_arrays_are_close(const Fiber::_2dArray<Matrix<ValueType> >& leftArrays,
                       const Fiber::_2dArray<Matrix<ValueType> >& rightArrays,
                       typename Fiber::ScalarTraits<ValueType>::RealType tolerance)
{
    typedef typename Fiber::ScalarTraits<ValueType>::RealType RealType;
    const int digits10 = std::numeric_limits<RealType>::digits10;

    boost::test_tools::predicate_result result(true);
    if (leftArrays.extent(0) != rightArrays.extent(0) ||
            leftArrays.extent(1) != rightArrays.extent(1)) {
        result = false;
        result.message() << "Size mismatch of Fiber::_2dArray [("
                         << leftArrays.extent(0) << ", "
                         << leftArrays.extent(1) << ") != ("
                         << rightArrays.extent(0) << ", "
                         << rightArrays.extent(1) << "]";
        return result;
    }

    for (size_t ra = 0; ra < leftArrays.extent(0); ++ra)
        for (size_t ca = 0; ca < leftArrays.extent(1); ++ca) {
            const Matrix<ValueType>& left = leftArrays(ra, ca);
            const Matrix<ValueType>& right = rightArrays(ra, ca);
            if (left.rows() != right.rows() || left.cols() != right.cols()) {
                result = false;
                result.message() << "Size mismatch of matrix ("
                                 << ra << ", " << ca << ") [("
                                 << left.rows() << ", " << left.cols() << ") != ("
                                 << right.rows() << ", " << right.cols() << ")]";
                return result;
            }
            for (size_t r = 0; r < left.rows(); ++r)
                for (size_t c = 0; c < left.cols(); ++c) {
                    RealType diff = std::abs(left(r, c) - right(r, c));
                    RealType avg = std::abs(left(r, c) + right(r, c)) / 2.;
                    if (diff > tolerance * (1 + avg)) {
                        result = false;
                        result.message() << std::setprecision(digits10)
                                         << "\n  Mismatch in matrix ("
                                         << ra << ", " << ca << ") at position ("
                                         << r << ", " << c << ") ["
                                         << left(r, c) << " != " << right(r, c) << "]";
                    }
                }
        }

    return result;
}



#endif
