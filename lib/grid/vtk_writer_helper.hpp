#ifndef bempp_vtk_writer_helper_hpp
#define bempp_vtk_writer_helper_hpp

#include <armadillo>
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
