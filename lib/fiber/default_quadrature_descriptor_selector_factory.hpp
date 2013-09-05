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

#ifndef fiber_default_quadrature_descriptor_selector_factory_hpp
#define fiber_default_quadrature_descriptor_selector_factory_hpp

#include "quadrature_descriptor_selector_factory.hpp"
#include "accuracy_options.hpp"

namespace Fiber
{

template <typename BasisFunctionType>
class DefaultQuadratureDescriptorSelectorFactory :
        public QuadratureDescriptorSelectorFactory<BasisFunctionType>
{
    typedef QuadratureDescriptorSelectorFactory<BasisFunctionType> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;

    explicit DefaultQuadratureDescriptorSelectorFactory(
        const AccuracyOptionsEx& accuracyOptions = AccuracyOptionsEx());

    virtual shared_ptr<QuadratureDescriptorSelectorForGridFunctions<CoordinateType> >
    makeQuadratureDescriptorSelectorForGridFunctions(
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType>*> >& testShapesets) const;

    virtual shared_ptr<QuadratureDescriptorSelectorForIntegralOperators<CoordinateType> >
    makeQuadratureDescriptorSelectorForIntegralOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType>*> >& trialShapesets) const;

    virtual shared_ptr<QuadratureDescriptorSelectorForLocalOperators<CoordinateType> >
    makeQuadratureDescriptorSelectorForLocalOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType>*> >& trialShapesets) const;

    virtual shared_ptr<QuadratureDescriptorSelectorForPotentialOperators<BasisFunctionType> >
    makeQuadratureDescriptorSelectorForPotentialOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType>*> >& trialShapesets) const;

private:
    AccuracyOptionsEx m_accuracyOptions;
};

} // namespace Fiber

#endif
