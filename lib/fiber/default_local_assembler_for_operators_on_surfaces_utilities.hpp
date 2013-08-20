// Copyright (C) 2011-2013 by the Bem++ Authors
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

#ifndef bempp_default_local_assembler_for_operators_utilities_hpp
#define bempp_default_local_assembler_for_operators_utilities_hpp

#include "../common/common.hpp"

#include "scalar_traits.hpp"
#include "../common/armadillo_fwd.hpp"

#include <vector>

namespace Fiber
{

template <typename BasisFunctionType> class Shapeset;
template <typename CoordinateType> class RawGridGeometry;

// These are internal routines used by the
// DefaultLocalAssemblerFor*OperatorsOnSurfaces classes. They are not part of
// the public interface of BEM++ and may change at any time.

template <typename BasisFunctionType>
struct DefaultLocalAssemblerForOperatorsOnSurfacesUtilities
{
    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    static void checkConsistencyOfGeometryAndShapesets(
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Shapeset<BasisFunctionType>*>& shapesets);

    static void precalculateElementSizesAndCentersForSingleGrid(
            const RawGridGeometry<CoordinateType>& rawGeometry,
            std::vector<CoordinateType>& elementSizesSquared,
            arma::Mat<CoordinateType>& elementCenters,
            CoordinateType& averageElementSize);

private:
    static CoordinateType elementSizeSquared(
            int elementIndex,
            const RawGridGeometry<CoordinateType>& rawGeometry);
    static arma::Col<CoordinateType> elementCenter(
            int elementIndex,
            const RawGridGeometry<CoordinateType>& rawGeometry);
};

} // namespace Fiber

#endif
