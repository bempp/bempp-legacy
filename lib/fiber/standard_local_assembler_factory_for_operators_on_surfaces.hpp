// Copyright (C) 2011-2012 by the Fiber Authors
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

#ifndef fiber_standard_local_assembler_factory_for_operators_on_surfaces_hpp
#define fiber_standard_local_assembler_factory_for_operators_on_surfaces_hpp

#include "local_assembler_factory.hpp"
#include "standard_local_assembler_for_identity_operator_on_surface.hpp"
#include "standard_local_assembler_for_integral_operators_on_surfaces.hpp"
#include "standard_local_assembler_for_source_terms_on_surfaces.hpp"
#include "opencl_options.hpp"

#include <stdexcept>

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
class StandardLocalAssemblerFactoryForOperatorsOnSurfaces :
        public LocalAssemblerFactory<ValueType, GeometryFactory>
{
private:
    typedef StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>
        LocalAssemblerForIntegralOperators_;
    typedef StandardLocalAssemblerForIdentityOperatorOnSurface<ValueType, GeometryFactory>
        LocalAssemblerForIdentityOperator_;
    typedef StandardLocalAssemblerForSourceTermsOnSurfaces<ValueType, GeometryFactory>
        LocalAssemblerForSourceTerms_;

public:
    virtual std::auto_ptr<LocalAssemblerForIntegralOperators<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& testExpression,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            const OpenClHandler& openClHandler,
            bool cacheSingularIntegrals) const {
        return std::auto_ptr<LocalAssemblerForIntegralOperators<ValueType> >(
                    new LocalAssemblerForIntegralOperators_(
                        geometryFactory, rawGeometry,
                        testBases, trialBases,
                        testExpression, kernel, trialExpression,
                        openClHandler, cacheSingularIntegrals));
    }

    virtual std::auto_ptr<LocalAssemblerForIntegralOperators<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            const OpenClHandler& openClHandler,
            bool cacheSingularIntegrals) const {
        throw std::runtime_error("StandardLocalAssemblerFactoryForOperatorsOnSurfaces::"
                                 "make(): collocation mode not implemented yet.");
    }

    virtual std::auto_ptr<LocalAssemblerForIdentityOperator<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& testExpression,
            const Expression<ValueType>& trialExpression,
            const OpenClHandler& openClHandler) const {
        return std::auto_ptr<LocalAssemblerForIdentityOperator<ValueType> >(
                    new LocalAssemblerForIdentityOperator_(
                        geometryFactory, rawGeometry,
                        testBases, trialBases,
                        testExpression, trialExpression,
                        openClHandler));
    }

    virtual std::auto_ptr<LocalAssemblerForIdentityOperator<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& trialExpression,
            const OpenClHandler& openClHandler) const {
    throw std::runtime_error("StandardLocalAssemblerFactoryForOperatorsOnSurfaces::"
                             "make(): collocation mode not implemented yet.");
    }

    virtual std::auto_ptr<LocalAssemblerForSourceTerms<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const Expression<ValueType>& testExpression,
            const Function<ValueType>& function,
            const OpenClHandler& openClHandler) const {
        return std::auto_ptr<LocalAssemblerForSourceTerms<ValueType> >(
                    new LocalAssemblerForSourceTerms_(
                        geometryFactory, rawGeometry,
                        testBases,
                        testExpression, function,
                        openClHandler));
    }
};

} // namespace Fiber

#endif
