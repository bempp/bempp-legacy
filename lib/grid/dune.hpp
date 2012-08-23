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

#ifndef bempp_dune_hpp
#define bempp_dune_hpp

#include "../common/common.hpp"

#include "bempp/common/config_alugrid.hpp"

#include <memory>
#include <stack> // fix a bug in foamgrid -- this header is not included where it should be
#include <dune/foamgrid/foamgrid.hh>
#ifdef WITH_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

namespace Bempp
{
typedef Dune::FoamGrid<3 /* dimWorld */> Default2dIn3dDuneGrid;
#ifdef WITH_ALUGRID
typedef Dune::ALUSimplexGrid<3 /*dimGrid */, 3 /* dimWorld */> Default3dIn3dDuneGrid;
#endif
} // namespace Bempp

#endif
