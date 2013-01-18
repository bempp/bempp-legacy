// Copyright (C) 2011-2013 by the BEM++ Authors
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

#ifndef bempp_precalculated_potential_operator_hpp
#define bempp_precalculated_potential_operator_hpp

#include "../common/common.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/shared_ptr.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename ValueType> class DiscreteBoundaryOperator;

template <typename BasisFunctionType, typename ResultType>
class PrecalculatedPotentialOperator
{
public:
    // name abuse: the discrete operator plays the part of a discrete *potential*
    // operator, not a discrete boundary operator...
    explicit PrecalculatedPotentialOperator(
            const shared_ptr<const DiscreteBoundaryOperator<ResultType> >& op);

    arma::Mat<ResultType> evaluate(
                const GridFunction<BasisFunctionType, ResultType>& argument) const;

private:
    /** \cond PRIVATE */
    const shared_ptr<const DiscreteBoundaryOperator<ResultType> > m_op;
    /** \endcond */
};

} // namespace Bempp

#endif
