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

#include "vtk_writer_helper.hpp"

#include "bempp/common/config_data_types.hpp"
#include "../common/scalar_traits.hpp"

#include <string>

namespace Bempp {

template <typename ResultType>
typename boost::enable_if<boost::is_complex<ResultType>, void>::type
exportSingleDataSetToVtk(VtkWriter &vtkWriter, const Matrix<ResultType> &data,
                         VtkWriter::DataType dataType, const char *dataLabel,
                         const char *fileNamesBase, const char *filesPath,
                         VtkWriter::OutputType outputType) {
  typedef typename ScalarTraits<ResultType>::RealType RealType;
  // TODO: figure out how to avoid this copy
  const Matrix<RealType> dataReal(data.real());
  const Matrix<RealType> dataImag(data.imag());
  const Matrix<RealType> dataAbs(data.cwiseAbs());

  if (dataType == VtkWriter::CELL_DATA) {
    vtkWriter.addCellData(dataReal, dataLabel + std::string(".r"));
    vtkWriter.addCellData(dataImag, dataLabel + std::string(".i"));
    vtkWriter.addCellData(dataAbs, dataLabel + std::string(".abs"));
  } else { // VERTEX_DATA
    vtkWriter.addVertexData(dataReal, dataLabel + std::string(".r"));
    vtkWriter.addVertexData(dataImag, dataLabel + std::string(".i"));
    vtkWriter.addVertexData(dataAbs, dataLabel + std::string(".abs"));
  }
  if (filesPath)
    vtkWriter.pwrite(fileNamesBase, filesPath, ".", outputType);
  else
    vtkWriter.write(fileNamesBase, outputType);
}

template <typename ResultType>
typename boost::disable_if<boost::is_complex<ResultType>, void>::type
exportSingleDataSetToVtk(VtkWriter &vtkWriter, const Matrix<ResultType> &data,
                         VtkWriter::DataType dataType, const char *dataLabel,
                         const char *fileNamesBase, const char *filesPath,
                         VtkWriter::OutputType outputType) {
  if (dataType == VtkWriter::CELL_DATA)
    vtkWriter.addCellData(data, dataLabel);
  else // VERTEX_DATA
    vtkWriter.addVertexData(data, dataLabel);
  if (filesPath)
    vtkWriter.pwrite(fileNamesBase, filesPath, ".", outputType);
  else
    vtkWriter.write(fileNamesBase, outputType);
}

#if defined(ENABLE_SINGLE_PRECISION)
template void
exportSingleDataSetToVtk(VtkWriter &vtkWriter, const Matrix<float> &data,
                         VtkWriter::DataType dataType, const char *dataLabel,
                         const char *fileNamesBase, const char *filesPath,
                         VtkWriter::OutputType outputType);

#if defined(ENABLE_COMPLEX_KERNELS) || defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
template void exportSingleDataSetToVtk(VtkWriter &vtkWriter,
                                       const Matrix<std::complex<float>> &data,
                                       VtkWriter::DataType dataType,
                                       const char *dataLabel,
                                       const char *fileNamesBase,
                                       const char *filesPath,
                                       VtkWriter::OutputType outputType);
#endif

#endif

#if defined(ENABLE_DOUBLE_PRECISION)
template void
exportSingleDataSetToVtk(VtkWriter &vtkWriter, const Matrix<double> &data,
                         VtkWriter::DataType dataType, const char *dataLabel,
                         const char *fileNamesBase, const char *filesPath,
                         VtkWriter::OutputType outputType);

#if defined(ENABLE_COMPLEX_KERNELS) || defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
template void exportSingleDataSetToVtk(VtkWriter &vtkWriter,
                                       const Matrix<std::complex<double>> &data,
                                       VtkWriter::DataType dataType,
                                       const char *dataLabel,
                                       const char *fileNamesBase,
                                       const char *filesPath,
                                       VtkWriter::OutputType outputType);
#endif

#endif

} // namespace Bempp
