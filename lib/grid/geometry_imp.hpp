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

#include "../common/common.hpp"

#include "geometry.hpp" // keep IDEs happy
#include "../fiber/geometrical_data.hpp"

#include <algorithm>

namespace Bempp {

inline void Geometry::setup(const Matrix<double> &corners,
                            const Vector<char> &auxData) {
  setupImpl(corners, auxData);
}

inline void Geometry::setup(const Matrix<float> &corners,
                            const Vector<char> &auxData) {
  Matrix<double> cornersDouble;
  convertMat(corners, cornersDouble);
  setupImpl(cornersDouble, auxData);
}

inline void Geometry::getCorners(Matrix<double> &c) const {
  getCornersImpl(c);
}

inline void Geometry::getCorners(Matrix<float> &c) const {
  Matrix<double> cDouble;
  getCornersImpl(cDouble);
  convertMat(cDouble, c);
}

inline void Geometry::local2global(const Matrix<double> &local,
                                   Matrix<double> &global) const {
  local2globalImpl(local, global);
}

inline void Geometry::local2global(const Matrix<float> &local,
                                   Matrix<float> &global) const {
  Matrix<double> localDouble, globalDouble;
  convertMat(local, localDouble);
  local2globalImpl(localDouble, globalDouble);
  convertMat(globalDouble, global);
}

inline void Geometry::global2local(const Matrix<double> &global,
                                   Matrix<double> &local) const {
  global2localImpl(global, local);
}

inline void Geometry::global2local(const Matrix<float> &global,
                                   Matrix<float> &local) const {
  Matrix<double> localDouble, globalDouble;
  convertMat(global, globalDouble);
  global2localImpl(globalDouble, localDouble);
  convertMat(localDouble, local);
}

inline void
Geometry::getIntegrationElements(const Matrix<double> &local,
                                 RowVector<double> &int_element) const {
  getIntegrationElementsImpl(local, int_element);
}

inline void
Geometry::getIntegrationElements(const Matrix<float> &local,
                                 RowVector<float> &int_element) const {
  Matrix<double> localDouble;
  RowVector<double> int_elementDouble;
  convertMat(local, localDouble);
  getIntegrationElementsImpl(localDouble, int_elementDouble);
  int_element.resize(int_elementDouble.cols());
  for (int i = 0; i < int_elementDouble.cols();++i)
      int_element(i) = int_elementDouble(i);
}

inline void Geometry::getCenter(Eigen::Ref<Vector<double>> c) const {
  getCenterImpl(c);
}

inline void Geometry::getCenter(Eigen::Ref<Vector<float>> c) const {

    Vector<double> cDouble;
    cDouble.resize(c.rows());
    getCenterImpl(cDouble);
    for (int i = 0; i < c.rows(); ++i)
        c(i) = cDouble(i);

}

inline void
Geometry::getJacobiansTransposed(const Matrix<double> &local,
                                 std::vector<Matrix<double>> &jacobian_t) const {
  const size_t n = local.cols();
  const size_t mdim = dim();
  const size_t cdim = dimWorld();
  jacobian_t.resize(n);
  getJacobiansTransposedImpl(local, jacobian_t);
}

inline void
Geometry::getJacobiansTransposed(const Matrix<float> &local,
                                 std::vector<Matrix<float>> &jacobian_t) const {
  Matrix<double> localDouble;
  std::vector<Matrix<double>> jacobian_tDouble;
  convertMat(local, localDouble);
  getJacobiansTransposedImpl(localDouble, jacobian_tDouble);
  jacobian_t.resize(jacobian_tDouble.size());
  for (int i = 0; i<jacobian_t.size();++i){
      jacobian_t[i].resize(jacobian_tDouble[i].rows(),jacobian_tDouble[i].cols());
      for (int j = 0; j < jacobian_tDouble[i].cols();++j)
          for (int k = 0; k < jacobian_tDouble[i].rows();++k)
              jacobian_t[i](k,j) = static_cast<float>(jacobian_tDouble[i](k,j));
  }

}

inline void
Geometry::getJacobiansTransposed(const Matrix<double> &local,
                                 Fiber::_3dArray<double> &jacobian_t) const {
  getJacobiansTransposedImpl(local, jacobian_t);
}

inline void
Geometry::getJacobiansTransposed(const Matrix<float> &local,
                                 Fiber::_3dArray<float> &jacobian_t) const {
  Matrix<double> localDouble;
  Fiber::_3dArray<double> jacobian_tDouble;
  convertMat(local, localDouble);
  getJacobiansTransposedImpl(localDouble, jacobian_tDouble);
  convertCube(jacobian_tDouble, jacobian_t);
}

inline void Geometry::getJacobianInversesTransposed(
    const Matrix<double> &local, std::vector<Matrix<double>> &jacobian_inv_t) const {
  getJacobianInversesTransposedImpl(local, jacobian_inv_t);
}

inline void Geometry::getJacobianInversesTransposed(
    const Matrix<float> &local, std::vector<Matrix<float>> &jacobian_inv_t) const {
  Matrix<double> localDouble;
  std::vector<Matrix<double>> jacobian_inv_tDouble;
  convertMat(local, localDouble);
  getJacobianInversesTransposed(localDouble, jacobian_inv_tDouble);
  jacobian_inv_t.resize(jacobian_inv_tDouble.size());
  for (int i = 0; i<jacobian_inv_t.size();++i){
      jacobian_inv_t[i].resize(jacobian_inv_tDouble[i].rows(),jacobian_inv_tDouble[i].cols());
      for (int j = 0; j < jacobian_inv_tDouble[i].cols();++j)
          for (int k = 0; k < jacobian_inv_tDouble[i].rows();++k)
              jacobian_inv_t[i](k,j) = static_cast<float>(jacobian_inv_tDouble[i](k,j));
  }

}

inline void Geometry::getJacobianInversesTransposed(
    const Matrix<double> &local,
    Fiber::_3dArray<double> &jacobian_inv_t) const {
  getJacobianInversesTransposedImpl(local, jacobian_inv_t);
}

inline void Geometry::getJacobianInversesTransposed(
    const Matrix<float> &local,
    Fiber::_3dArray<float> &jacobian_inv_t) const {
  Matrix<double> localDouble;
  Fiber::_3dArray<double> jacobian_inv_tDouble;
  convertMat(local, localDouble);
  getJacobianInversesTransposedImpl(localDouble, jacobian_inv_tDouble);
  convertCube(jacobian_inv_tDouble, jacobian_inv_t);
}

inline void Geometry::getNormals(const Matrix<double> &local,
                                 Matrix<double> &normal) const {
  getNormalsImpl(local, normal);
}

inline void Geometry::getNormals(const Matrix<float> &local,
                                 Matrix<float> &normal) const {
  Matrix<double> localDouble, normalDouble;
  convertMat(local, localDouble);
  getNormalsImpl(localDouble, normalDouble);
  convertMat(normalDouble, normal);
}

inline void Geometry::getData(size_t what, const Matrix<double> &local,
                              Fiber::GeometricalData<double> &data) const {
  getDataImpl(what, local, data);
}

inline void Geometry::getData(size_t what, const Matrix<float> &local,
                              Fiber::GeometricalData<float> &data) const {
  Matrix<double> localDouble;
  Fiber::GeometricalData<double> dataDouble;
  convertMat(local, localDouble);
  getDataImpl(what, localDouble, dataDouble);
  convertMat(dataDouble.globals, data.globals);
  convertMat(dataDouble.normals, data.normals);
  convertMat(dataDouble.integrationElements, data.integrationElements);
  convertCube(dataDouble.jacobiansTransposed, data.jacobiansTransposed);
  convertCube(dataDouble.jacobianInversesTransposed,
              data.jacobianInversesTransposed);
}

template <typename T1, typename T2>
void Geometry::convertMat(const Matrix<T1> &in, Matrix<T2> &out) const {
  out.resize(in.rows(), in.cols());
  for (int j = 0; j < in.cols();++j )
      for (int i = 0 ; i < in.rows(); ++i)
          out(i,j) = in(i,j);
}

template <typename T1, typename T2>
void Geometry::convertMat(const Vector<T1> &in, Vector<T2> &out) const {
  out.resize(in.rows());
  for (int j = 0; j < in.rows();++j )
      out(j) = in(j);
}

template <typename T1, typename T2>
void Geometry::convertMat(const RowVector<T1> &in, RowVector<T2> &out) const {
  out.resize(in.cols());
  for (int j = 0; j < in.cols();++j )
      out(j) = in(j);
}



template <typename T1, typename T2>
void Geometry::convertCube(const Fiber::_3dArray<T1> &in,
                           Fiber::_3dArray<T2> &out) const {
  out.set_size(in.extent(0), in.extent(1), in.extent(2));
  std::copy(in.begin(), in.end(), out.begin());
}

} // namespace Bempp

#endif
