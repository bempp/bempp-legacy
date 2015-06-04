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
// and has been relicensed with permission by the HyENA authors. This does not
// affect the license of any
// other part of HyENA.

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
template <typename T> struct ScalarTraits {
  enum { is_scalar = false };
};

template <> struct ScalarTraits<int> {
  enum { is_scalar = true };
};

template <> struct ScalarTraits<float> {
  enum { is_scalar = true };
};

template <> struct ScalarTraits<double> {
  enum { is_scalar = true };
};

template <> struct ScalarTraits<std::complex<float>> {
  enum { is_scalar = true };
};

template <> struct ScalarTraits<std::complex<double>> {
  enum { is_scalar = true };
};

#endif // include guard
