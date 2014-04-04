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

#ifndef bempp_test_entity_hpp
#define bempp_test_entity_hpp

#include "test_grid.hpp"
#include "grid/entity.hpp"
#include "grid/entity_pointer.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/grid_view.hpp"

struct TriangularEntityManager : public SimpleTriangularGridManager {
    TriangularEntityManager() : SimpleTriangularGridManager() {
    }

    template <int codim>
    typename std::auto_ptr<Bempp::EntityPointer<codim> > getPointerToSecondEntityOnLevel0() {
        std::auto_ptr<Bempp::GridView> bemppGridView = bemppGrid->levelView(0);
        std::auto_ptr<Bempp::EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
        it->next();
        return std::auto_ptr<Bempp::EntityPointer<codim> >(it.release());
    }

    template <int codim>
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator getDunePointerToSecondEntityOnLevel0() {
        DuneGrid::LevelGridView duneGridView = duneGrid->levelView(0);
        typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneIt = duneGridView.begin<codim>();
        ++duneIt;
        return duneIt;
    }
};

#endif
