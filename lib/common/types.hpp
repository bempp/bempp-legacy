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

#include "../fiber/types.hpp"
#include "../grid/common.hpp" // to be merged with this file
#include "../grid/geometry_type.hpp"
#include "../grid/index_set.hpp"

namespace Bempp
{

using Fiber::CallVariant;
using Fiber::TEST_TRIAL;
using Fiber::TRIAL_TEST;

using Fiber::ALL_DOFS;

typedef int ElementVariant;

struct EntityIndex
{
    GeometryType type;
    IndexSet::IndexType index;

    EntityIndex() {
    }

    EntityIndex(const GeometryType& type_, IndexSet::IndexType index_) :
        type(type_), index(index_) {
    }

    bool operator==(const EntityIndex& other) const {
        return (type == other.type && index == other.index);
    }

    bool operator<(const EntityIndex& other) const {
        if (type == other.type)
            return index < other.index;
        else
            return type < other.type;
    }
};

typedef int GlobalDofIndex;
typedef Fiber::LocalDofIndex LocalDofIndex;

struct LocalDof
{
    LocalDof() {}
    LocalDof(EntityIndex ei, LocalDofIndex ldi) :
        entityIndex(ei), dofIndex(ldi) {
    }

    EntityIndex entityIndex;
    LocalDofIndex dofIndex;
};

struct Point3D
{
    ctype x, y, z;
};

} // namespace Bempp

#endif
