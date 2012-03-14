// Copyright (C) 2011-2012 by the Fiber Authors
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

// File to be removed once unit tests are in place

#include "standard_local_assembler_factory_for_operators_on_surfaces.hpp"
#include "../grid/geometry.hpp"
#include "../grid/geometry_factory.hpp"

namespace Fiber
{

#ifdef COMPILE_FOR_FLOAT
template class StandardLocalAssemblerFactoryForOperatorsOnSurfaces<float, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class StandardLocalAssemblerFactoryForOperatorsOnSurfaces<std::complex<float>, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class StandardLocalAssemblerFactoryForOperatorsOnSurfaces<std::complex<double>, Bempp::GeometryFactory>;
#endif

} // namespace Fiber
