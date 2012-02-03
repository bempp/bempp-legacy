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

#ifndef bempp_quadrature_selector_factory_hpp
#define bempp_quadrature_selector_factory_hpp

#include <memory>

namespace Bempp
{

class AssemblyOptions;
template <typename ValueType> class QuadratureSelector;

template <typename ValueType>
class QuadratureSelectorFactory
{
public:
    static std::auto_ptr<QuadratureSelector<ValueType> > make(
            const AssemblyOptions& options);
};

/*
To choose quadrature rule, we need:
- polynomial orders of elements in all directions
- positions of elements -- in practice, their vertices
- whether elements are curved or not, and if yes, polynomial order of that curvature
- kernel properties (regular/singular) and extra quadrature order necessary
  because of the kernel's variation (note: this is dependent on the elements'
  size, since the kernel's variation follows a scale of its own!)
*/

} //namespace Bempp

#endif
