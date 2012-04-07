// Copyright (C) 2011 by the BEM++ Authors
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

#ifndef bempp_meshes_hpp
#define bempp_meshes_hpp

#include <memory>

namespace Bempp
{
    class Grid;
} // namespace Bempp

enum MeshVariant
{
    TWO_DISJOINT_TRIANGLES,
    TWO_TRIANGLES_SHARING_VERTEX_0,
    TWO_TRIANGLES_SHARING_VERTICES_2_AND_0,
    TWO_TRIANGLES_SHARING_VERTICES_1_AND_0,
    TWO_TRIANGLES_SHARING_EDGES_0_AND_0,
    TWO_TRIANGLES_SHARING_EDGES_1_AND_0,
    SIMPLE_MESH_9,
    CUBE_12,
    CUBE_12_REORIENTED, // all elements oriented so that normals point outwards
    CUBE_384,
    CUBE_6144,
    CUBE_24576,
    SPHERE_152,
    SPHERE_644
};

std::auto_ptr<Bempp::Grid> loadMesh(MeshVariant mv);

void dumpElementList(const Bempp::Grid* grid);

#endif
