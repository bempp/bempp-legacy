#ifndef fiber_geometrical_data_hpp
#define fiber_geometrical_data_hpp

#include <armadillo>

namespace Fiber
{

enum GeometricalDataType
{
    GLOBALS = 0x0001,
    INTEGRATION_ELEMENTS = 0x0002,
    NORMALS = 0x0004,
    JACOBIANS_TRANSPOSED = 0x0008,
    JACOBIAN_INVERSES_TRANSPOSED = 0x0010
};

template <typename CoordinateType>
struct GeometricalData
{
    arma::Mat<CoordinateType> globals;
    arma::Row<CoordinateType> integrationElements;
    arma::Cube<CoordinateType> jacobiansTransposed;
    arma::Cube<CoordinateType> jacobianInversesTransposed;
};

} // namespace Fiber

#endif // GEOMETRY_DATA_TYPES_HPP
