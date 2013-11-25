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

#include "modified_helmholtz_3d_calderon_projector.hpp"
#include "../common/shared_ptr.hpp"

#include "modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
#include "identity_operator.hpp"
#include "blocked_operator_structure.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <boost/make_shared.hpp>


namespace Bempp {

template <typename BasisFunctionType, typename KernelType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
modifiedHelmholtz3dExteriorCalderonProjector(
        const shared_ptr<const Context<BasisFunctionType,ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& hminusSpace,
        const shared_ptr<const Space<BasisFunctionType> >& hplusSpace,
        KernelType waveNumber,
        const std::string& label,
        int symmetry,
        bool useInterpolation,
        int interpPtsPerWavelength) {

    typedef BoundaryOperator<BasisFunctionType,ResultType> BdOp;

    shared_ptr<const Space<BasisFunctionType> > internalHplusSpace = hplusSpace->discontinuousSpace(hplusSpace);

    BdOp internalSlp = modifiedHelmholtz3dSingleLayerBoundaryOperator(context,internalHplusSpace,internalHplusSpace,internalHplusSpace,waveNumber,label+"_slp",symmetry,
                                                                          useInterpolation,interpPtsPerWavelength);

    BdOp dlp = modifiedHelmholtz3dDoubleLayerBoundaryOperator(context,hplusSpace,hplusSpace,hminusSpace,waveNumber,label+"_dlp",symmetry,
                                                                                  useInterpolation,interpPtsPerWavelength);

    BdOp adjDlp = adjoint(dlp);

    BdOp hyp = modifiedHelmholtz3dHypersingularBoundaryOperator(context,hplusSpace,hminusSpace,hplusSpace,waveNumber,label+"_hyp",symmetry,
                                                                        useInterpolation,interpPtsPerWavelength,internalSlp);

    BdOp idSpaceTransformation1 = identityOperator(context,hminusSpace,internalHplusSpace,internalHplusSpace);
    BdOp idSpaceTransformation2 = identityOperator(context,internalHplusSpace,hplusSpace,hminusSpace);
    BdOp idDouble = identityOperator(context,hplusSpace,hplusSpace,hminusSpace);
    BdOp idAdjDouble = identityOperator(context,hminusSpace,hminusSpace,hplusSpace);

    // Now Assemble the entries of the Calderon Projector

    BlockedOperatorStructure<BasisFunctionType,ResultType> structure;

    structure.setBlock(0,0,.5*idDouble+dlp);
    structure.setBlock(0,1,-1.*idSpaceTransformation2*internalSlp*idSpaceTransformation1);
    structure.setBlock(1,0,-1.*hyp);
    structure.setBlock(1,1,.5*idAdjDouble-adjDlp);

    return BlockedBoundaryOperator<BasisFunctionType,ResultType>(structure);



}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS, KERNEL, RESULT) \
    template BlockedBoundaryOperator<BASIS, RESULT> \
    modifiedHelmholtz3dExteriorCalderonProjector( \
        const shared_ptr<const Context<BASIS, RESULT> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        KERNEL, \
        const std::string&, int, bool, int)

FIBER_ITERATE_OVER_BASIS_KERNEL_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);


} // Bempp
