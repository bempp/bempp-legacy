// Copyright (C) 2011-2012 by the BEM++ Authors
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

#ifndef fiber_single_quadrature_rule_family_hpp
#define fiber_single_quadrature_rule_family_hpp

#include "../common/common.hpp"
#include "../common/armadillo_fwd.hpp"
#include "single_quadrature_descriptor.hpp"

#include <vector>

namespace Fiber {

/** \brief Family of quadrature rules over single elements. */
template <typename CoordinateType> class SingleQuadratureRuleFamily {
public:
  /** \brief Destructor. */
  virtual ~SingleQuadratureRuleFamily() {}

  /** \brief Fill arrays of quadrature points and weights.
   *
   *  \param[in]  desc    Quadrature descriptor.
   *  \param[out] points  2D array whose (i, j)th element will
   *    contain the ith (local) coordinate of the jth quadrature point.
   *  \param[out] weights Vector of quadrature point weights. */
  virtual void fillQuadraturePointsAndWeights(
      const SingleQuadratureDescriptor &desc, arma::Mat<CoordinateType> &points,
      std::vector<CoordinateType> &weights) const = 0;
};

} // namespace Fiber

#endif
