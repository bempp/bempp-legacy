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


#include "interpolated_function.hpp"

#include "../fiber/geometrical_data.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/vtk_writer.hpp"
#include "../grid/vtk_writer_helper.hpp"

#include "../space/piecewise_linear_continuous_scalar_space.hpp"

namespace Bempp
{

template <typename ValueType>
InterpolatedFunction<ValueType>::InterpolatedFunction(
        const Grid& grid, const arma::Mat<ValueType>& vertexValues,
        InterpolationMethod method) :
    m_grid(grid), m_vertexValues(vertexValues), m_method(method)
{
    std::unique_ptr<GridView> view = grid.leafView();

    if (view->entityCount(grid.dim()) != vertexValues.n_cols)
        throw std::invalid_argument("VolumeGridFunction::VolumeGridFunction(): "
                                    "dimension of expansionCoefficients "
                                    "does not match the number of grid vertices");
    if (method != LINEAR)
        throw std::invalid_argument("VolumeGridFunction::VolumeGridFunction(): "
                                    "invalid interpolation method");
}

template <typename ValueType>
const Grid& InterpolatedFunction<ValueType>::grid() const
{
    return m_grid;
}

template <typename ValueType>
int InterpolatedFunction<ValueType>::worldDimension() const
{
    return m_grid.dimWorld();
}

template <typename ValueType>
int InterpolatedFunction<ValueType>::codomainDimension() const
{
    return m_vertexValues.n_rows;
}

template <typename ValueType>
void InterpolatedFunction<ValueType>::addGeometricalDependencies(size_t& geomDeps) const
{
    geomDeps |= Fiber::GLOBALS;
}

template <typename ValueType>
void InterpolatedFunction<ValueType>::evaluate(
        const Fiber::GeometricalData<CoordinateType>& geomData,
        arma::Mat<ValueType>& result) const
{

#ifndef NDEBUG
    const arma::Mat<CoordinateType>& points = geomData.globals;
    if ((int)points.n_rows != worldDimension())
        throw std::invalid_argument("InterpolatedFunction::evaluate(): "
                                    "incompatible world dimension");
#endif

    throw std::runtime_error("InterpolatedFunction::evaluate(): "
                             "not implemented yet");
}

template <typename ValueType>
void InterpolatedFunction<ValueType>::exportToVtk(
        const char* dataLabel, const char* fileNamesBase, const char* filesPath,
        VtkWriter::OutputType outputType) const
{
    std::unique_ptr<GridView> view = m_grid.leafView();
    std::unique_ptr<VtkWriter> vtkWriter = view->vtkWriter();

    exportSingleDataSetToVtk(*vtkWriter, m_vertexValues, VtkWriter::VERTEX_DATA,
                             dataLabel, fileNamesBase, filesPath, outputType);
}

//template <typename ValueType>
//void InterpolatedFunction<ValueType>::setSurfaceValues(
//        const GridFunction<ValueType>& surfaceFunction)
//{
//    throw std::runtime_error("InterpolatedFunction::setSurfaceValues(): "
//                             "not implemented yet");
//}

//template <typename ValueType>
//void InterpolatedFunction<ValueType>::setSurfaceValues(
//        const InterpolatedFunction<ValueType>& surfaceFunction)
//{
//    throw std::runtime_error("InterpolatedFunction::setSurfaceValues(): "
//                             "not implemented yet");
//}

template <typename ValueType>
void InterpolatedFunction<ValueType>::checkCompatibility(
        const InterpolatedFunction<ValueType> &other) const
{
    if (&m_grid != &other.m_grid ||
            m_vertexValues.n_rows != other.m_vertexValues.n_rows ||
            m_vertexValues.n_cols != other.m_vertexValues.n_cols ||
            m_method != other.m_method)
        throw std::runtime_error("InterpolatedFunction::checkCompatibility(): "
                                 "incompatible operands");
}

template <typename ValueType>
InterpolatedFunction<ValueType>&
InterpolatedFunction<ValueType>::operator+=(
        const InterpolatedFunction<ValueType> &rhs)
{
    checkCompatibility(rhs);
    m_vertexValues += rhs.m_vertexValues;
    return *this;
}

template <typename ValueType>
InterpolatedFunction<ValueType>&
InterpolatedFunction<ValueType>::operator-=(
        const InterpolatedFunction<ValueType> &rhs)
{
    checkCompatibility(rhs);
    m_vertexValues -= rhs.m_vertexValues;
    return *this;
}

template <typename ValueType>
InterpolatedFunction<ValueType>&
InterpolatedFunction<ValueType>::operator*=(ValueType rhs)
{
    m_vertexValues *= rhs;
    return *this;
}

template <typename ValueType>
InterpolatedFunction<ValueType>&
InterpolatedFunction<ValueType>::operator/=(ValueType rhs)
{
    m_vertexValues *= static_cast<ValueType>(1.) / rhs;
    return *this;
}

template <typename ValueType>
const InterpolatedFunction<ValueType>
InterpolatedFunction<ValueType>::operator+(
        const InterpolatedFunction<ValueType> &other) const
{
    return InterpolatedFunction<ValueType>(*this) += other;
}

template <typename ValueType>
const InterpolatedFunction<ValueType>
InterpolatedFunction<ValueType>::operator-(
        const InterpolatedFunction<ValueType> &other) const
{
    return InterpolatedFunction<ValueType>(*this) -= other;
}

template <typename ValueType>
const InterpolatedFunction<ValueType>
InterpolatedFunction<ValueType>::operator/(ValueType other) const
{
    return InterpolatedFunction<ValueType>(*this) /= other;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(InterpolatedFunction);

} // namespace Bempp
