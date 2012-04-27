#include "vtk_writer_helper.hpp"

#include "config_data_types.hpp"
#include "../common/scalar_traits.hpp"

#include <string>

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
        VtkWriter::OutputType outputType)
{
    typedef typename ScalarTraits<ResultType>::RealType RealType;
    // TODO: figure out how to avoid this copy
    const arma::Mat<RealType> dataReal(arma::real(data));
    const arma::Mat<RealType> dataImag(arma::imag(data));

    if (dataType == VtkWriter::CELL_DATA) {
        vtkWriter.addCellData(dataReal, dataLabel + std::string(".r"));
        vtkWriter.addCellData(dataImag, dataLabel + std::string(".i"));
    } else { // VERTEX_DATA
        vtkWriter.addVertexData(dataReal, dataLabel + std::string(".r"));
        vtkWriter.addVertexData(dataImag, dataLabel + std::string(".i"));
    }
    if (filesPath)
        vtkWriter.pwrite(fileNamesBase, filesPath, ".", outputType);
    else
        vtkWriter.write(fileNamesBase, outputType);
}

template <typename ResultType>
typename boost::disable_if<boost::is_complex<ResultType>, void>::type
exportSingleDataSetToVtk(
        VtkWriter& vtkWriter,
        const arma::Mat<ResultType>& data,
        VtkWriter::DataType dataType,
        const char* dataLabel,
        const char* fileNamesBase, const char* filesPath,
        VtkWriter::OutputType outputType)
{
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
template void exportSingleDataSetToVtk(
VtkWriter& vtkWriter,
const arma::Mat<float>& data,
VtkWriter::DataType dataType,
const char* dataLabel,
const char* fileNamesBase, const char* filesPath,
VtkWriter::OutputType outputType);

#  if defined(ENABLE_COMPLEX_KERNELS) || defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
template void exportSingleDataSetToVtk(
VtkWriter& vtkWriter,
const arma::Mat<std::complex<float> >& data,
VtkWriter::DataType dataType,
const char* dataLabel,
const char* fileNamesBase, const char* filesPath,
VtkWriter::OutputType outputType);
#  endif

#endif

#if defined(ENABLE_DOUBLE_PRECISION)
template void exportSingleDataSetToVtk(
VtkWriter& vtkWriter,
const arma::Mat<double>& data,
VtkWriter::DataType dataType,
const char* dataLabel,
const char* fileNamesBase, const char* filesPath,
VtkWriter::OutputType outputType);

#  if defined(ENABLE_COMPLEX_KERNELS) || defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
template void exportSingleDataSetToVtk(
VtkWriter& vtkWriter,
const arma::Mat<std::complex<double> >& data,
VtkWriter::DataType dataType,
const char* dataLabel,
const char* fileNamesBase, const char* filesPath,
VtkWriter::OutputType outputType);
#  endif

#endif

} // namespace Bempp
