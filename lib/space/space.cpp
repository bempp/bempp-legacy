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

#include "space.hpp"
#include "../common/config_data_types.hpp"

#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"

#include <complex>

namespace Bempp
{

/** \brief Get pointers to Basis objects corresponding to all elements. */
template <typename BasisFunctionType>
void getAllBases(const Space<BasisFunctionType>& space,
        std::vector<const Fiber::Basis<BasisFunctionType>*>& bases)
{
    std::auto_ptr<GridView> view = space.grid().leafView();
    const int elementCount = view->entityCount(0);

    bases.clear();
    bases.reserve(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished()) {
        bases.push_back(&space.basis(it->entity()));
        it->next();
    }
}

#define INSTANTIATE(TYPE) \
    template \
    void getAllBases(const Space<TYPE>& space, \
                std::vector<const Fiber::Basis<TYPE>*>& bases)

#ifdef ENABLE_SINGLE_PRECISION
INSTANTIATE(float);
#  ifdef ENABLE_COMPLEX_BASIS_FUNCTIONS
INSTANTIATE(std::complex<float>);
#  endif
#endif

#ifdef ENABLE_DOUBLE_PRECISION
INSTANTIATE(double);
#  ifdef ENABLE_COMPLEX_BASIS_FUNCTIONS
INSTANTIATE(std::complex<double>);
#  endif
#endif

} // namespace Bempp
