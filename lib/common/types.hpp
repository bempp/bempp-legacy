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

#ifndef bempp_types_hpp
#define bempp_types_hpp

/** \file common/types.hpp Types and headers used in various parts of BEM++ */
#include "../fiber/types.hpp"
#include "../grid/geometry_type.hpp"
#include "../grid/index_set.hpp"
#include "common.hpp"

namespace Bempp {

using Fiber::CallVariant;
using Fiber::TEST_TRIAL;
using Fiber::TRIAL_TEST;

using Fiber::ALL_DOFS;

typedef boost::property_tree::ptree ParameterList;

typedef int ElementVariant;
typedef int EntityIndex;
typedef int GlobalDofIndex;
typedef int FlatLocalDofIndex;
typedef Fiber::LocalDofIndex LocalDofIndex;

/** \ingroup weak_form_assembly_internal
 *  \brief Local degree of freedom. */
struct LocalDof {
  LocalDof() {}
  LocalDof(EntityIndex ei, LocalDofIndex ldi)
      : entityIndex(ei), dofIndex(ldi) {}

  EntityIndex entityIndex;
  LocalDofIndex dofIndex;
};

/** \ingroup weak_form_assembly_internal
 *  \brief Point in a three-dimensional space. */
template <typename ValueType> struct Point3D { ValueType x, y, z; };
template <typename ValueType>
ValueType &accessPointByIndex(Point3D<ValueType> &point, int i) {
  if (i == 0)
    return point.x;
  if (i == 1)
    return point.y;
  if (i == 2)
    return point.z;
}

template <typename ValueType>
const ValueType &accessPointByIndex(const Point3D<ValueType> &point, int i) {
  if (i == 0)
    return point.x;
  if (i == 1)
    return point.y;
  if (i == 2)
    return point.z;
}

} // namespace Bempp

#endif
