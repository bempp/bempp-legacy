#ifndef bempp_local_assembler_construction_helper_hpp
#define bempp_local_assembler_construction_helper_hpp

#include "assembly_options.hpp"
#include "../common/not_implemented_error.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../space/space.hpp"

#include <boost/make_shared.hpp>

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
