#include "numerical_double_integrator.hpp" // To keep IDEs happy

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
void NumericalDoubleIntegrator<ValueType, GeometryFactory>::setupGeometry(
        int elementIndex,
        const arma::Mat<ValueType>& vertices,
        const arma::Mat<int>& elementCornerIndices,
        const arma::Mat<char>& auxElementData,
        typename GeometryFactory::Geometry& geometry)
{
    const int dimGrid = vertices.n_rows;
    int cornerCount = 0;
    for (; cornerCount < elementCornerIndices.n_rows; ++cornerCount)
        if (elementCornerIndices(cornerCount, elementIndex) < 0)
            break;
    arma::Mat<ValueType> corners(dimGrid, cornerCount);
    for (int cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex)
        corners.col(cornerIndex) = vertices.col(
                    elementCornerIndices(cornerIndex, elementIndex));
    geometry.setup(corners, auxElementData.unsafe_col(elementIndex));
}

} // namespace Fiber
