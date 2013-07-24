#ifndef bempp_vc_geometrical_data_hpp
#define bempp_vc_geometrical_data_hpp

#include <Vc/Vc>

namespace Fiber
{

template <typename CoordinateType>
struct VcGeometricalData
{
    VcGeometricalData(int pointCount_, int dimWorld) :
        pointCount(pointCount_),
        paddedPointCount(((pointCount + Vc::Vector<CoordinateType>::Size - 1) / 
                          Vc::Vector<CoordinateType>::Size) * Vc::Vector<CoordinateType>::Size),
        globals(paddedPointCount * dimWorld)
        {}

    size_t pointCount;
    size_t paddedPointCount;
    // One-dimensional array of vectors. Possibly could be replaced by a Fiber::_1dArray
    Vc::Memory<Vc::Vector<CoordinateType>, 0u, 0u> globals;
};

}

#endif
