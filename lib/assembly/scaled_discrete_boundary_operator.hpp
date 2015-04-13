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

#include "bempp/common/config_ahmed.hpp"
#include "bempp/common/config_trilinos.hpp"

#ifndef bempp_scaled_discrete_boundary_operator_hpp
#define bempp_scaled_discrete_boundary_operator_hpp

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include "discrete_boundary_operator.hpp"

#include "../common/shared_ptr.hpp"

#ifdef WITH_TRILINOS
#include <Teuchos_RCP.hpp>
#endif

namespace Bempp {

/** \ingroup composite_discrete_boundary_operators
 *  \brief Scaled discrete boundary operator.
 *
 *  This class represents a discrete boundary operator multiplied by a scalar.
 */
template <typename ValueType>
class ScaledDiscreteBoundaryOperator
    : public DiscreteBoundaryOperator<ValueType> {
public:
  typedef DiscreteBoundaryOperator<ValueType> Base;

  /** \brief Constructor.
   *
   *  Construct the discrete boundary operator \f$\alpha L\f$, where
   *  \f$\alpha\f$ is the scalar \p multiplier and \f$L\f$ is the operator
   *  represented by \p op. */
  ScaledDiscreteBoundaryOperator(ValueType multiplier,
                                 const shared_ptr<const Base> &op);

  virtual Matrix<ValueType> asMatrix() const;

  virtual unsigned int rowCount() const;
  virtual unsigned int columnCount() const;

  virtual void addBlock(const std::vector<int> &rows,
                        const std::vector<int> &cols, const ValueType alpha,
                        Matrix<ValueType> &block) const;

#ifdef WITH_AHMED
  shared_ptr<const DiscreteBoundaryOperator<ValueType>>
  asDiscreteAcaBoundaryOperator(double eps = -1, int maximumRank = -1,
                                bool interleave = false) const;
#endif

private:
  virtual void applyBuiltInImpl(const TranspositionMode trans,
                                const Eigen::Ref<Vector<ValueType>> &x_in,
                                Eigen::Ref<Vector<ValueType>> y_inout,
                                const ValueType alpha,
                                const ValueType beta) const;

private:
  ValueType m_multiplier;
  shared_ptr<const Base> m_operator;
};

} // namespace Bempp

#endif
