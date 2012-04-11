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

#ifndef bempp_test_index_set_hpp
#define bempp_test_index_set_hpp

#include "test_grid_view.hpp"

struct TriangularLevel0IndexSetManager : public TriangularLevel0GridViewManager {
    TriangularLevel0IndexSetManager() : TriangularLevel0GridViewManager(),
        bemppIndexSet(bemppGridView->indexSet()),
        duneIndexSet(duneGridView.indexSet()) {
    }

    typedef Bempp::Default2dIn3dDuneGrid::LevelGridView::IndexSet DuneIndexSet;
    const Bempp::IndexSet& bemppIndexSet;
    const DuneIndexSet& duneIndexSet;
};

struct TriangularLeafIndexSetManager : public TriangularLeafGridViewManager {
    TriangularLeafIndexSetManager() : TriangularLeafGridViewManager(),
        bemppIndexSet(bemppGridView->indexSet()),
        duneIndexSet(duneGridView.indexSet()) {
    }

    typedef Bempp::Default2dIn3dDuneGrid::LeafGridView::IndexSet DuneIndexSet;
    const Bempp::IndexSet& bemppIndexSet;
    const DuneIndexSet& duneIndexSet;
};

#endif
