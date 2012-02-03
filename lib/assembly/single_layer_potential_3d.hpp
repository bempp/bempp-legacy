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

#ifndef bempp_single_layer_potential_3d_hpp
#define bempp_single_layer_potential_3d_hpp

#include "elementary_weakly_singular_integral_operator.hpp"
#include "single_layer_potential_3d_kernel.hpp"
#include "../space/all_shape_functions.hpp"

namespace Bempp
{

template <typename ValueType>
class SingleLayerPotential3D :
        public ElementaryWeaklySingularIntegralOperator<ValueType>
{
private:
    virtual const Kernel<ValueType>& kernel() const {
        return m_kernel;
    }

    virtual std::auto_ptr<FunctionFamily<ValueType> > testFunctionFamily(
            const Space<ValueType>& testSpace,
            ElementVariant elementVariant) const {
        return anyFunctionFamily(testSpace, elementVariant);
    }

    virtual std::auto_ptr<FunctionFamily<ValueType> > trialFunctionFamily(
            const Space<ValueType>& trialSpace,
            ElementVariant elementVariant) const {
        return anyFunctionFamily(trialSpace, elementVariant);
    }

private:
    virtual std::auto_ptr<FunctionFamily<ValueType> > anyFunctionFamily(
            const Space<ValueType>& space,
            ElementVariant elementVariant) const {
        typedef AllShapeFunctions<ValueType> Family;
        return std::auto_ptr<FunctionFamily<ValueType> >(
                    new Family(space, elementVariant));
    }

    SingleLayerPotential3DKernel<ValueType> m_kernel;
};

} // namespace Bempp

#endif
