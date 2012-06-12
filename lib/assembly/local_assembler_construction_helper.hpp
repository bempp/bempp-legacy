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

#ifndef bempp_local_assembler_construction_helper_hpp
#define bempp_local_assembler_construction_helper_hpp

#include "../common/common.hpp"

#include "assembly_options.hpp"
#include "../common/not_implemented_error.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../space/space.hpp"

#include "../common/boost_make_shared_fwd.hpp"

namespace Bempp
{

struct LocalAssemblerConstructionHelper
{
    template <typename CoordinateType>
    static void collectGridData(
            const Grid& grid,
            shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& rawGeometry,
            shared_ptr<GeometryFactory>& geometryFactory) {
        typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;

        rawGeometry = boost::make_shared<RawGridGeometry>(grid.dim(),
                                                          grid.dimWorld());
        std::auto_ptr<GridView> view = grid.leafView();
        view->getRawElementData(
                    rawGeometry->vertices(), rawGeometry->elementCornerIndices(),
                    rawGeometry->auxData());
        geometryFactory = shared_ptr<GeometryFactory>(
                grid.elementGeometryFactory().release());
    }

    template <typename BasisFunctionType>
    static void collectBases(
            const Space<BasisFunctionType>& space,
            shared_ptr<std::vector<const Fiber::Basis<BasisFunctionType>*> >& bases) {
        typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;

        if (!space.dofsAssigned())
            throw std::runtime_error(
                    "LocalAssemblerConstructionHelper::collectBases(): "
                    "degrees of freedom must be assigned "
                    "before calling this function");
        bases = boost::make_shared<BasisPtrVector>();
        getAllBases(space, *bases);
    }

    // Probably in future will be generalised to arbitrary number of grids

    template <typename CoordinateType>
    static void makeOpenClHandler(
            const OpenClOptions& openClOptions,
            const shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& rawGeometry,
            shared_ptr<Fiber::OpenClHandler>& openClHandler) {
        openClHandler = boost::make_shared<Fiber::OpenClHandler>(openClOptions);
        if (openClHandler->UseOpenCl())
            openClHandler->pushGeometry(rawGeometry->vertices(),
                                        rawGeometry->elementCornerIndices());
    }

    template <typename CoordinateType>
    static void makeOpenClHandler(
            const OpenClOptions& openClOptions,
            const shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            shared_ptr<Fiber::OpenClHandler>& openClHandler) {
        openClHandler = boost::make_shared<Fiber::OpenClHandler>(openClOptions);
        if (openClHandler->UseOpenCl()) {
            if (testRawGeometry.get() == trialRawGeometry.get())
                openClHandler->pushGeometry(testRawGeometry->vertices(),
                                            testRawGeometry->elementCornerIndices());
            else
                throw NotImplementedError(
                        "LocalAssemblerConstructionHelper::makeOpenClHandler(): "
                        "OpenCL-based assembly of operators with test and trial "
                        "spaces defined on different grids is not currently supported");
        }
    }
};

} // namespace Bempp

#endif
