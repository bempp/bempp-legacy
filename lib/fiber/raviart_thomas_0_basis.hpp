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

#ifndef fiber_raviart_thomas_0_basis_hpp
#define fiber_raviart_thomas_0_basis_hpp

#include "../common/common.hpp"

#include "basis.hpp"

#include "basis_data.hpp"
#include "dune_basis_helper.hpp"

//#include <dune/localfunctions/raviartthomas.hh> ///raviartthomas0q2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas02d/raviartthomas02dlocalbasis.hh>

namespace Fiber
{

template <int elementVertexCount, typename CoordinateType, typename ValueType>
struct RaviartThomas0BasisTraits
{
};

// Triangle
template <typename CoordinateType, typename ValueType>
struct RaviartThomas0BasisTraits<3, CoordinateType, ValueType>
{
public:
    typedef Dune::RT02DLocalBasis<CoordinateType, ValueType> DuneBasis;
};

// // Quadrilateral
// template <typename CoordinateType, typename ValueType>
// struct RaviartThomasOrder0BasisTraits<4, CoordinateType, ValueType>
// {
// public:
//     typedef Dune::Q1LocalBasis<CoordinateType, ValueType, 2> DuneBasis;
// };

template <int elementVertexCount, typename ValueType>
class RaviartThomas0Basis : public Basis<ValueType>
{
public:
    typedef typename Basis<ValueType>::CoordinateType CoordinateType;

private:
    typedef typename RaviartThomas0BasisTraits
    <elementVertexCount, CoordinateType, ValueType>::DuneBasis DuneBasis;

public:
    virtual int size() const {
        DuneBasis basis;
        return basis.size();
    }

    virtual int order() const {
        return 1;
    }

    virtual void evaluate(size_t what,
                          const arma::Mat<CoordinateType>& points,
                          LocalDofIndex localDofIndex,
                          BasisData<ValueType>& data) const {
        if (localDofIndex != ALL_DOFS &&
                (localDofIndex < 0 || size() <= localDofIndex))
            throw std::invalid_argument("RaviartThomas0Basis::"
                                        "evaluate(): Invalid localDofIndex");

        if (what & VALUES)
            evaluateBasisFunctionsWithDune<CoordinateType, ValueType, DuneBasis>(
                        points, localDofIndex, data.values);
        if (what & DERIVATIVES)
            evaluateBasisFunctionDerivativesWithDune<CoordinateType, ValueType, DuneBasis>(
                        points, localDofIndex, data.derivatives);
    }
};

} // namespace Fiber

#endif
