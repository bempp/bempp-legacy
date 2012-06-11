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
#include "fiber/array_2d.hpp"
#include "fiber/array_4d.hpp"

#include <armadillo>
#include <complex>
#include <boost/test/unit_test.hpp>

template <typename ValueType>
boost::test_tools::predicate_result
check_arrays_are_close(const arma::Mat<ValueType>& left,
                       const arma::Mat<ValueType>& right,
                       typename Fiber::ScalarTraits<ValueType>::RealType tolerance)
{
    boost::test_tools::predicate_result result(true);
    if (left.n_rows != right.n_rows || left.n_cols != right.n_cols) {
        result = false;
        result.message() << "Size mismatch [("
                         << left.n_rows << ", " << left.n_cols << ") != ("
                         << right.n_rows << ", " << right.n_cols << ")]";
        return result;
    }
    for (size_t r = 0; r < left.n_rows; ++r)
        for (size_t c = 0; c < left.n_cols; ++c) {
            typename Fiber::ScalarTraits<ValueType>::RealType diff =
                    std::abs(left(r, c) - right(r, c));
            typename Fiber::ScalarTraits<ValueType>::RealType avg =
                    std::abs(left(r, c) + right(r, c)) / 2.;
            if (diff > tolerance * (1 + avg)) {
                result = false;
                result.message() << "\n  Mismatch at position ("
                                 << r << ", " << c << ") ["
                                 << left(r, c) << " != " << right(r, c) << "]";
            }
        }

    return result;
}

template <typename ValueType>
boost::test_tools::predicate_result
check_arrays_are_close(const arma::Cube<ValueType>& left,
                       const arma::Cube<ValueType>& right,
                       typename Fiber::ScalarTraits<ValueType>::RealType tolerance)
{
    boost::test_tools::predicate_result result(true);
    if (left.n_rows != right.n_rows ||
            left.n_cols != right.n_cols ||
            left.n_slices != right.n_slices) {
        result = false;
        result.message() << "Size mismatch [("
                         << left.n_rows << ", " << left.n_cols << ", "
                         << left.n_slices << ") != ("
                         << right.n_rows << ", " << right.n_cols << ", "
                         << right.n_slices << ")]";
        return result;
    }
    for (size_t r = 0; r < left.n_rows; ++r)
        for (size_t c = 0; c < left.n_cols; ++c)
            for (size_t s = 0; s < left.n_slices; ++s) {
                typename Fiber::ScalarTraits<ValueType>::RealType diff =
                        std::abs(left(r, c, s) - right(r, c, s));
                typename Fiber::ScalarTraits<ValueType>::RealType avg =
                        std::abs(left(r, c, s) + right(r, c, s)) / 2.;
                if (diff > tolerance * (1 + avg)) {
                    result = false;
                    result.message() << "\n  Mismatch at position ("
                                     << r << ", " << c << ", " << s << ") ["
                                     << left(r, c, s) << " != "
                                     << right(r, c, s) << "]";
                }
            }

    return result;
}

template <typename ValueType>
boost::test_tools::predicate_result
check_arrays_are_close(const Fiber::Array4d<ValueType>& left,
                       const Fiber::Array4d<ValueType>& right,
                       typename Fiber::ScalarTraits<ValueType>::RealType tolerance)
{
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
                    typename Fiber::ScalarTraits<ValueType>::RealType diff =
                            std::abs(left(r, c, s, u) - right(r, c, s, u));
                    typename Fiber::ScalarTraits<ValueType>::RealType avg =
                            std::abs(left(r, c, s, u) + right(r, c, s, u)) / 2.;
                    if (diff > tolerance * (1 + avg)) {
                        result = false;
                        result.message() << "\n  Mismatch at position ("
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
check_arrays_are_close(const Fiber::Array2d<arma::Mat<ValueType> >& leftArrays,
                       const Fiber::Array2d<arma::Mat<ValueType> >& rightArrays,
                       typename Fiber::ScalarTraits<ValueType>::RealType tolerance)
{
    boost::test_tools::predicate_result result(true);
    if (leftArrays.extent(0) != rightArrays.extent(0) ||
            leftArrays.extent(1) != rightArrays.extent(1)) {
        result = false;
        result.message() << "Size mismatch of Fiber::Array2d [("
                         << leftArrays.extent(0) << ", "
                         << leftArrays.extent(1) << ") != ("
                         << rightArrays.extent(0) << ", "
                         << rightArrays.extent(1) << "]";
        return result;
    }

    for (size_t ra = 0; ra < leftArrays.extent(0); ++ra)
        for (size_t ca = 0; ca < leftArrays.extent(1); ++ca) {
            const arma::Mat<ValueType>& left = leftArrays(ra, ca);
            const arma::Mat<ValueType>& right = rightArrays(ra, ca);
            if (left.n_rows != right.n_rows || left.n_cols != right.n_cols) {
                result = false;
                result.message() << "Size mismatch of matrix ("
                                 << ra << ", " << ca << ") [("
                                 << left.n_rows << ", " << left.n_cols << ") != ("
                                 << right.n_rows << ", " << right.n_cols << ")]";
                return result;
            }
            for (size_t r = 0; r < left.n_rows; ++r)
                for (size_t c = 0; c < left.n_cols; ++c) {
                    typename Fiber::ScalarTraits<ValueType>::RealType diff =
                            std::abs(left(r, c) - right(r, c));
                    typename Fiber::ScalarTraits<ValueType>::RealType avg =
                            std::abs(left(r, c) + right(r, c)) / 2.;
                    if (diff > tolerance * (1 + avg)) {
                        result = false;
                        result.message() << "\n  Mismatch in matrix ("
                                         << ra << ", " << ca << ") at position ("
                                         << r << ", " << c << ") ["
                                         << left(r, c) << " != " << right(r, c) << "]";
                    }
                }
        }

    return result;
}

template <typename ValueType>
boost::test_tools::predicate_result
check_arrays_are_close(const std::vector<arma::Mat<ValueType> >& leftArrays,
                       const std::vector<arma::Mat<ValueType> >& rightArrays,
                       typename Fiber::ScalarTraits<ValueType>::RealType tolerance)
{
    boost::test_tools::predicate_result result(true);
    if (leftArrays.size() != rightArrays.size()) {
        result = false;
        result.message() << "Size mismatch of std::vector ["
                         << leftArrays.size() << " != "
                         << rightArrays.size() << "]";
        return result;
    }

    for (size_t ra = 0; ra < leftArrays.size(); ++ra) {
        const arma::Mat<ValueType>& left = leftArrays[ra];
        const arma::Mat<ValueType>& right = rightArrays[ra];
        if (left.n_rows != right.n_rows || left.n_cols != right.n_cols) {
            result = false;
            result.message() << "Size mismatch of matrix ("
                             << ra << ") [("
                             << left.n_rows << ", " << left.n_cols << ") != ("
                             << right.n_rows << ", " << right.n_cols << ")]";
            return result;
        }
        for (size_t r = 0; r < left.n_rows; ++r)
            for (size_t c = 0; c < left.n_cols; ++c) {
                typename Fiber::ScalarTraits<ValueType>::RealType diff =
                        std::abs(left(r, c) - right(r, c));
                typename Fiber::ScalarTraits<ValueType>::RealType avg =
                        std::abs(left(r, c) + right(r, c)) / 2.;
                if (diff > tolerance * (1 + avg)) {
                    result = false;
                    result.message() << "\n  Mismatch in matrix ("
                                     << ra << ") at position ("
                                     << r << ", " << c << ") ["
                                     << left(r, c) << " != " << right(r, c) << "]";
                }
            }
    }

    return result;
}

#endif
