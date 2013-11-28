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

#ifndef bempp_fmm_modified_helmholtz_3d_high_frequency_hpp
#define bempp_fmm_modified_helmholtz_3d_high_frequency_hpp

#include "fmm_transform.hpp"

#include <vector>

#include "../common/armadillo_fwd.hpp"
#include "../fiber/scalar_traits.hpp"

// must include interpolate_on_sphere since cannot have vectors of incomplete types
#include "interpolate_on_sphere.hpp"

namespace Bempp
{

/** \cond FORWARD_DECL */
//template <typename ValueType> class InterpolateOnSphere;
//template <typename ValueType> class FmmTransform;
/** \endcond */

// all operations are diagonal in the case of plane wave expansion
template <typename ValueType>
class FmmModifiedHelmholtz3dHighFrequency : public FmmTransform<ValueType>
{
public:
    typedef typename FmmTransform<ValueType>::CoordinateType CoordinateType;

    FmmModifiedHelmholtz3dHighFrequency(ValueType kappa, unsigned int expansionOrder, 
        unsigned int expansionOrderMax, unsigned int levels);

    // multipole to multipole (M2M) translation matrix
    virtual arma::Mat<ValueType> M2M(
        const arma::Col<CoordinateType>& childPosition, 
        const arma::Col<CoordinateType>& parentPosition,
        unsigned int level) const;

    // multipole to local (M2L) translation matrix
    virtual arma::Mat<ValueType> M2L(
        const arma::Col<CoordinateType>& sourceCentre, 
        const arma::Col<CoordinateType>& fieldCentre,
        const arma::Col<CoordinateType>& boxSize,
        unsigned int level) const;

    // local to local (L2L) translation matrix
    virtual arma::Mat<ValueType> L2L(
        const arma::Col<CoordinateType>& parentPosition, 
        const arma::Col<CoordinateType>& childPosition,
        unsigned int level) const;

    // interpolate multipole and local coefficients between levels
    virtual void interpolate(
        unsigned int levelOld,
        unsigned int levelNew,
        const arma::Col<ValueType>& coefficientsOld, 
        arma::Col<ValueType>& coefficientsNew) const;

    virtual void generateGaussPoints();

    //virtual void getKernelWeight(arma::Mat<ValueType>& kernelWeight) const;

    ValueType kappa() const {return m_kappa;}
private:
    ValueType m_kappa;
    std::vector<unsigned int> m_Ls;

    std::vector<InterpolateOnSphere<ValueType> > m_interpolatorsUpwards;
    std::vector<InterpolateOnSphere<ValueType> > m_interpolatorsDownwards;
};

} // namespace Bempp

#endif
