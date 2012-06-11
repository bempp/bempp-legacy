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

#ifndef bempp_geometry_imp_hpp
#define bempp_geometry_imp_hpp

#include "geometry.hpp" // keep IDEs happy
#include "../fiber/geometrical_data.hpp"

namespace Bempp
{

inline void Geometry::setup(const arma::Mat<double>& corners,
                            const arma::Col<char>& auxData)
{
    setupImpl(corners, auxData);
}

inline void Geometry::setup(const arma::Mat<float>& corners,
                            const arma::Col<char>& auxData)
{
    arma::Mat<double> cornersDouble;
    convertMat(corners, cornersDouble);
    setupImpl(cornersDouble, auxData);
}

inline void Geometry::getCorners(arma::Mat<double>& c) const
{
    getCornersImpl(c);
}

inline void Geometry::getCorners(arma::Mat<float>& c) const
{
    arma::Mat<double> cDouble;
    getCornersImpl(cDouble);
    convertMat(cDouble, c);
}

inline void Geometry::local2global(const arma::Mat<double>& local,
                                   arma::Mat<double>& global) const
{
    local2globalImpl(local, global);
}

inline void Geometry::local2global(const arma::Mat<float>& local,
                                   arma::Mat<float>& global) const
{
    arma::Mat<double> localDouble, globalDouble;
    convertMat(local, localDouble);
    local2globalImpl(localDouble, globalDouble);
    convertMat(globalDouble, global);
}

inline void Geometry::global2local(const arma::Mat<double>& global,
                                   arma::Mat<double>& local) const
{
    global2localImpl(global, local);
}

inline void Geometry::global2local(const arma::Mat<float>& global,
                                   arma::Mat<float>& local) const
{
    arma::Mat<double> localDouble, globalDouble;
    convertMat(global, globalDouble);
    global2localImpl(globalDouble, localDouble);
    convertMat(localDouble, local);
}

inline void Geometry::getIntegrationElements(const arma::Mat<double>& local,
                                             arma::Row<double>& int_element) const
{
    getIntegrationElementsImpl(local, int_element);
}

inline void Geometry::getIntegrationElements(const arma::Mat<float>& local,
                                             arma::Row<float>& int_element) const
{
    arma::Mat<double> localDouble;
    arma::Row<double> int_elementDouble;
    convertMat(local, localDouble);
    getIntegrationElementsImpl(localDouble, int_elementDouble);
    convertMat(int_elementDouble, int_element);
}

inline void Geometry::getCenter(arma::Col<double>& c) const
{
    getCenterImpl(c);
}

inline void Geometry::getCenter(arma::Col<float>& c) const
{
    arma::Col<double> cDouble;
    getCenterImpl(cDouble);
    convertMat(cDouble, c);
}

inline void Geometry::getJacobiansTransposed(
        const arma::Mat<double>& local,
        arma::Cube<double>& jacobian_t) const
{
    getJacobiansTransposedImpl(local, jacobian_t);
}

inline void Geometry::getJacobiansTransposed(
        const arma::Mat<float>& local,
        arma::Cube<float>& jacobian_t) const
{
    arma::Mat<double> localDouble;
    arma::Cube<double> jacobian_tDouble;
    convertMat(local, localDouble);
    getJacobiansTransposedImpl(localDouble, jacobian_tDouble);
    convertCube(jacobian_tDouble, jacobian_t);
}

inline void Geometry::getJacobianInversesTransposed(
        const arma::Mat<double>& local,
        arma::Cube<double>& jacobian_inv_t) const
{
    getJacobianInversesTransposedImpl(local, jacobian_inv_t);
}

inline void Geometry::getJacobianInversesTransposed(
        const arma::Mat<float>& local,
        arma::Cube<float>& jacobian_inv_t) const
{
    arma::Mat<double> localDouble;
    arma::Cube<double> jacobian_inv_tDouble;
    convertMat(local, localDouble);
    getJacobianInversesTransposedImpl(localDouble, jacobian_inv_tDouble);
    convertCube(jacobian_inv_tDouble, jacobian_inv_t);
}

inline void Geometry::getNormals(const arma::Mat<double>& local,
                                 arma::Mat<double>& normal) const
{
    getNormalsImpl(local, normal);
}

inline void Geometry::getNormals(const arma::Mat<float>& local,
                                 arma::Mat<float>& normal) const
{
    arma::Mat<double> localDouble, normalDouble;
    convertMat(local, localDouble);
    getNormalsImpl(localDouble, normalDouble);
    convertMat(normalDouble, normal);
}

inline void Geometry::getData(size_t what, const arma::Mat<double>& local,
                              Fiber::GeometricalData<double>& data) const
{
    getDataImpl(what, local, data);
}

inline void Geometry::getData(size_t what, const arma::Mat<float>& local,
                              Fiber::GeometricalData<float>& data) const
{
    arma::Mat<double> localDouble;
    Fiber::GeometricalData<double> dataDouble;
    convertMat(local, localDouble);
    getDataImpl(what, localDouble, dataDouble);
    convertMat(dataDouble.globals, data.globals);
    convertMat(dataDouble.normals, data.normals);
    convertMat(dataDouble.integrationElements, data.integrationElements);
    convertCube(dataDouble.jacobiansTransposed, data.jacobiansTransposed);
    convertCube(dataDouble.jacobiansTransposed, data.jacobiansTransposed);
}

template <typename T1, typename T2>
void Geometry::convertMat(const arma::Mat<T1>& in, arma::Mat<T2>& out) const
{
    out.set_size(in.n_rows, in.n_cols);
    for (size_t elem = 0; elem < in.n_elem; ++elem)
        out[elem] = in[elem];
}

template <typename T1, typename T2>
void Geometry::convertCube(const arma::Cube<T1>& in, arma::Cube<T2>& out) const
{
    out.set_size(in.n_rows, in.n_cols, in.n_slices);
    for (size_t elem = 0; elem < in.n_elem; ++elem)
        out[elem] = in[elem];
}

} // namespace Bempp

#endif
