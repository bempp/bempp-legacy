// Copyright (C) 2011-2015 by the BEM++ Authors
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

#include "hmat_interface.hpp"

#include "../hmat/block_cluster_tree.hpp"
#include "../hmat/geometry_data_type.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../common/types.hpp"

namespace Bempp {

template <typename BasisFunctionType>
SpaceHMatGeometryInterface<BasisFunctionType>::SpaceHMatGeometryInterface(
    const Space<BasisFunctionType> &space)
    : m_counter(0) {
  space.getGlobalDofBoundingBoxes(m_bemppBoundingBoxes);
}

template <typename BasisFunctionType>
shared_ptr<const hmat::GeometryDataType>
SpaceHMatGeometryInterface<BasisFunctionType>::next() {

  if (m_counter == m_bemppBoundingBoxes.size())
    return shared_ptr<hmat::GeometryDataType>();

  auto lbound = m_bemppBoundingBoxes[m_counter].lbound;
  auto ubound = m_bemppBoundingBoxes[m_counter].ubound;
  auto center = m_bemppBoundingBoxes[m_counter].reference;
  m_counter++;
  return shared_ptr<hmat::GeometryDataType>(new hmat::GeometryDataType(
      hmat::BoundingBox(lbound.x, ubound.x, lbound.y, ubound.y, lbound.z,
                        ubound.z),
      std::array<double, 3>({{center.x, center.y, center.z}})));
}

template <typename BasisFunctionType>
std::size_t
SpaceHMatGeometryInterface<BasisFunctionType>::numberOfEntities() const {
  return m_bemppBoundingBoxes.size();
}

template <typename BasisFunctionType>
void SpaceHMatGeometryInterface<BasisFunctionType>::reset() {
  m_counter = 0;
}

template <typename BasisFunctionType>
shared_ptr<hmat::DefaultBlockClusterTreeType>
generateBlockClusterTree(const Space<BasisFunctionType> &testSpace,
                         const Space<BasisFunctionType> &trialSpace,
                         const ParameterList &parameterList) {

  hmat::Geometry testGeometry;
  hmat::Geometry trialGeometry;

  auto testSpaceGeometryInterface = shared_ptr<hmat::GeometryInterface>(
      new SpaceHMatGeometryInterface<BasisFunctionType>(testSpace));

  auto trialSpaceGeometryInterface = shared_ptr<hmat::GeometryInterface>(
      new SpaceHMatGeometryInterface<BasisFunctionType>(trialSpace));

  hmat::fillGeometry(testGeometry, *testSpaceGeometryInterface);
  hmat::fillGeometry(trialGeometry, *trialSpaceGeometryInterface);

  auto admissibility =
      parameterList.template get<std::string>("options.hmat.admissibility");
  auto minBlockSize =
      parameterList.template get<int>("options.hmat.minBlockSize");
  auto maxBlockSize =
      parameterList.template get<int>("options.hmat.maxBlockSize");

  hmat::AdmissibilityFunction admissibilityFunction;

  if (admissibility == "strong") {
    auto eta = parameterList.template get<double>("options.hmat.eta");
    admissibilityFunction = hmat::StrongAdmissibility(eta);
  } else if (admissibility == "weak") {
    admissibilityFunction = hmat::WeakAdmissibility();
  } else if (admissibility == "high_freq") {
    auto waveNumber = parameterList.template get<double>("options.hmat.waveNumber");
    admissibilityFunction = hmat::HighFrequencyAdmissibility(waveNumber);
  } else
    throw std::runtime_error(
        "generateBlockClusterTree(): Unknown admissibility type");

  auto testClusterTree = shared_ptr<hmat::DefaultClusterTreeType>(
      new hmat::DefaultClusterTreeType(testGeometry, minBlockSize));

  auto trialClusterTree = shared_ptr<hmat::DefaultClusterTreeType>(
      new hmat::DefaultClusterTreeType(trialGeometry, minBlockSize));

  shared_ptr<hmat::DefaultBlockClusterTreeType> blockClusterTree(
      new hmat::DefaultBlockClusterTreeType(testClusterTree, trialClusterTree,
                                            maxBlockSize,
                                            admissibilityFunction));

  return blockClusterTree;
}

#define INSTANTIATE_NONMEMBER_FUNCTION(VALUE)                                  \
  template shared_ptr<hmat::DefaultBlockClusterTreeType>                       \
  generateBlockClusterTree(const Space<VALUE> &testSpace,                      \
                           const Space<VALUE> &trialSpace,                     \
                           const ParameterList &parameterList);

FIBER_ITERATE_OVER_VALUE_TYPES(INSTANTIATE_NONMEMBER_FUNCTION);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(SpaceHMatGeometryInterface);
}
