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
#ifndef modified_helmholtz_3d_calderon_projector_hpp
#define modified_helmholtz_3d_calderon_projector_hpp


#include "blocked_boundary_operator.hpp"
#include "../fiber/scalar_traits.hpp"
#include "helmholtz_3d_operators_common.hpp"
#include "symmetry.hpp"



namespace Bempp {

template <typename BasisFunctionType, typename ResultType>  class Context;
template <typename BasisFunctionType> class Space;

template <typename BasisFunctionType, typename KernelType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
modifiedHelmholtz3dExteriorCalderonProjector(
        const shared_ptr<const Context<BasisFunctionType,ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& hminusSpace,
        const shared_ptr<const Space<BasisFunctionType> >& hplusSpace,
        KernelType waveNumber,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY,
        bool useInterpolation = false,
        int interpPtsPerWavelength = DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY);

template <typename BasisFunctionType, typename KernelType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
modifiedHelmholtz3dInteriorCalderonProjector(
        const shared_ptr<const Context<BasisFunctionType,ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& hminusSpace,
        const shared_ptr<const Space<BasisFunctionType> >& hplusSpace,
        KernelType waveNumber,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY,
        bool useInterpolation = false,
        int interpPtsPerWavelength = DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY);


} // Bempp


#endif
