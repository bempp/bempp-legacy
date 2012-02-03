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

#include "../grid/geometry_type.hpp"
#include "../grid/common.hpp" // to be merged with this file

namespace Bempp
{

// typedef int ElementVariant;

struct ElementVariant
{
    unsigned char order[2];
    unsigned char vertexCount;
    unsigned char reserved;

    bool operator==(const ElementVariant& other) const {
        return (order[0] == other.order[0] &&
                order[1] == other.order[1] &&
                vertexCount == other.vertexCount &&
                reserved == other.reserved);
    }
    bool operator<(const ElementVariant& other) const {
        return (order[0] < other.order[0] ? true :
                (order[1] < other.order[1] ? true :
                (vertexCount < other.vertexCount ? true :
                (reserved < other.reserved))));
    }
};

struct EntityIndex
{
    GeometryType type;
    int index;

    EntityIndex(const GeometryType& type_, int index_) :
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
typedef int LocalDofIndex;

struct LocalDof
{
    EntityIndex entityIndex;
    LocalDofIndex dofIndex;
};

struct Point3D
{
    ctype x, y, z;
};

} // namespace Bempp

#endif
