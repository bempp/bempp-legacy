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

#ifndef bempp_vtk_writer_helper_hpp
#define bempp_vtk_writer_helper_hpp

#include "../common/common.hpp"

#include "../common/armadillo_fwd.hpp"
#include <boost/type_traits/is_complex.hpp>
#include <boost/utility/enable_if.hpp>

#include "vtk_writer.hpp"

namespace Bempp
{

template <typename ResultType>
typename boost::enable_if<boost::is_complex<ResultType>, void>::type
exportSingleDataSetToVtk(
        VtkWriter& vtkWriter,
        const arma::Mat<ResultType>& data,
        VtkWriter::DataType dataType,
        const char* dataLabel,
        const char* fileNamesBase, const char* filesPath,
        VtkWriter::OutputType outputType);

template <typename ResultType>
typename boost::disable_if<boost::is_complex<ResultType>, void>::type
exportSingleDataSetToVtk(
        VtkWriter& vtkWriter,
        const arma::Mat<ResultType>& data,
        VtkWriter::DataType dataType,
        const char* dataLabel,
        const char* fileNamesBase, const char* filesPath,
        VtkWriter::OutputType outputType);

} // namespace Bempp

#endif
