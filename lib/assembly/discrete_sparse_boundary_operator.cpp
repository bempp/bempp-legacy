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

#include "bempp/common/config_trilinos.hpp"
#ifdef WITH_TRILINOS

#include "discrete_sparse_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>
#include <stdexcept>

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialComm.h>
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>

namespace Bempp
{

// Helper functions for the applyBuiltIn member function
namespace
{

template <typename ValueType>
void reallyApplyBuiltInImpl(const Epetra_CrsMatrix& mat,
                            const TranspositionMode trans,
                            const arma::Col<ValueType>& x_in,
                            arma::Col<ValueType>& y_inout,
                            const ValueType alpha,
                            const ValueType beta);

template <>
void reallyApplyBuiltInImpl<double>(const Epetra_CrsMatrix& mat,
                                    const TranspositionMode trans,
                                    const arma::Col<double>& x_in,
                                    arma::Col<double>& y_inout,
                                    const double alpha,
                                    const double beta)
{
    if (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE) {
        assert(mat.NumGlobalRows() == static_cast<int>(x_in.n_rows));
        assert(mat.NumGlobalCols() == static_cast<int>(y_inout.n_rows));
    } else {
        assert(mat.NumGlobalCols() == static_cast<int>(x_in.n_rows));
        assert(mat.NumGlobalRows() == static_cast<int>(y_inout.n_rows));
    }

    Epetra_Map map_x(x_in.n_rows, 0, Epetra_SerialComm());
    Epetra_Map map_y(y_inout.n_rows, 0, Epetra_SerialComm());

    Epetra_Vector vec_x(View, map_x, const_cast<double*>(x_in.memptr()));
    // vec_temp will store the result of matrix * x_in
    Epetra_Vector vec_temp(map_y, false /* no need to initialise to zero */);

    mat.Multiply(trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE,
                 vec_x, vec_temp);

    if (beta == 0.)
        for (size_t i = 0; i < y_inout.n_rows; ++i)
            y_inout(i) = alpha * vec_temp[i];
    else
        for (size_t i = 0; i < y_inout.n_rows; ++i)
            y_inout(i) = alpha * vec_temp[i] + beta * y_inout(i);
}

template <>
void reallyApplyBuiltInImpl<float>(const Epetra_CrsMatrix& mat,
                                   const TranspositionMode trans,
                                   const arma::Col<float>& x_in,
                                   arma::Col<float>& y_inout,
                                   const float alpha,
                                   const float beta)
{
    // Copy the float vectors to double vectors
    arma::Col<double> x_in_double(x_in.n_rows);
    std::copy(x_in.begin(), x_in.end(), x_in_double.begin());
    arma::Col<double> y_inout_double(y_inout.n_rows);
    if (beta != 0.f)
        std::copy(y_inout.begin(), y_inout.end(), y_inout_double.begin());

    // Do the operation on the double vectors
    reallyApplyBuiltInImpl<double>(
                mat, trans, x_in_double, y_inout_double, alpha, beta);

    // Copy the result back to the float vector
    std::copy(y_inout_double.begin(), y_inout_double.end(), y_inout.begin());
}

template <>
void reallyApplyBuiltInImpl<std::complex<float> >(
        const Epetra_CrsMatrix& mat,
        const TranspositionMode trans,
        const arma::Col<std::complex<float> >& x_in,
        arma::Col<std::complex<float> >& y_inout,
        const std::complex<float> alpha,
        const std::complex<float> beta)
{
    // Do the y_inout *= beta part
    const std::complex<float> zero(0.f, 0.f);
    if (beta == zero) 
        y_inout.fill(zero);
    else
        y_inout *= beta;

    // Separate the real and imaginary components and store them in
    // double-precision vectors
    arma::Col<double> x_real(x_in.n_rows);
    for (size_t i = 0; i < x_in.n_rows; ++i)
        x_real(i) = x_in(i).real();
    arma::Col<double> x_imag(x_in.n_rows);
    for (size_t i = 0; i < x_in.n_rows; ++i)
        x_imag(i) = x_in(i).imag();
    arma::Col<double> y_real(y_inout.n_rows);
    for (size_t i = 0; i < y_inout.n_rows; ++i)
        y_real(i) = y_inout(i).real();
    arma::Col<double> y_imag(y_inout.n_rows);
    for (size_t i = 0; i < y_inout.n_rows; ++i)
        y_imag(i) = y_inout(i).imag();

    // Do the "+= alpha A x" part (in steps)
    reallyApplyBuiltInImpl<double>(
                mat, trans, x_real, y_real, alpha.real(), 1.);
    reallyApplyBuiltInImpl<double>(
                mat, trans, x_imag, y_real, -alpha.imag(), 1.);
    reallyApplyBuiltInImpl<double>(
                mat, trans, x_real, y_imag, alpha.imag(), 1.);
    reallyApplyBuiltInImpl<double>(
                mat, trans, x_imag, y_imag, alpha.real(), 1.);

    // Copy the result back to the complex vector
    for (size_t i = 0; i < y_inout.n_rows; ++i)
        y_inout(i) = std::complex<float>(y_real(i), y_imag(i));
}

template <>
void reallyApplyBuiltInImpl<std::complex<double> >(
        const Epetra_CrsMatrix& mat,
        const TranspositionMode trans,
        const arma::Col<std::complex<double> >& x_in,
        arma::Col<std::complex<double> >& y_inout,
        const std::complex<double> alpha,
        const std::complex<double> beta)
{
    // Do the y_inout *= beta part
    const std::complex<double> zero(0., 0.);
    if (beta == zero) 
        y_inout.fill(zero);
    else
        y_inout *= beta;

    // Separate the real and imaginary components
    arma::Col<double> x_real(arma::real(x_in));
    arma::Col<double> x_imag(arma::imag(x_in));
    arma::Col<double> y_real(arma::real(y_inout));
    arma::Col<double> y_imag(arma::imag(y_inout));

    // Do the "+= alpha A x" part (in steps)
    reallyApplyBuiltInImpl<double>(
                mat, trans, x_real, y_real, alpha.real(), 1.);
    reallyApplyBuiltInImpl<double>(
                mat, trans, x_imag, y_real, -alpha.imag(), 1.);
    reallyApplyBuiltInImpl<double>(
                mat, trans, x_real, y_imag, alpha.imag(), 1.);
    reallyApplyBuiltInImpl<double>(
                mat, trans, x_imag, y_imag, alpha.real(), 1.);

    // Copy the result back to the complex vector
    for (size_t i = 0; i < y_inout.n_rows; ++i)
        y_inout(i) = std::complex<double>(y_real(i), y_imag(i));
}

} // namespace


template <typename ValueType>
DiscreteSparseBoundaryOperator<ValueType>::
DiscreteSparseBoundaryOperator(const shared_ptr<const Epetra_CrsMatrix>& mat,
                             Symmetry symmetry, TranspositionMode trans) :
    m_mat(mat), m_symmetry(symmetry), m_trans(trans)
{
    m_domainSpace = Thyra::defaultSpmdVectorSpace<ValueType>(
                isTransposed() ? m_mat->NumGlobalRows() : m_mat->NumGlobalCols());
    m_rangeSpace = Thyra::defaultSpmdVectorSpace<ValueType>(
                isTransposed() ? m_mat->NumGlobalCols() : m_mat->NumGlobalRows());
}

template <typename ValueType>
void DiscreteSparseBoundaryOperator<ValueType>::dump() const
{
    if (isTransposed())
        std::cout << "Transpose of " << *m_mat << std::endl;
    else
        std::cout << *m_mat << std::endl;
}

template <typename ValueType>
arma::Mat<ValueType>
DiscreteSparseBoundaryOperator<ValueType>::asMatrix() const
{
    if (m_mat->Comm().NumProc() != 1)
        throw std::runtime_error(
                "DiscreteSparseBoundaryOperator::asMatrix(): "
                "conversion of distributed matrices to local matrices is unsupported");

    bool transposed = isTransposed();
    const int untransposedRowCount = m_mat->NumGlobalRows();
    arma::Mat<ValueType> mat(rowCount(), columnCount());
    mat.fill(0.);
    for (int row = 0; row < untransposedRowCount; ++row)
    {
        int entryCount = 0;
        double* values = 0;
        int* indices = 0;
        int errorCode = m_mat->ExtractMyRowView(row, entryCount,
                                                values, indices);
        if (errorCode != 0)
            throw std::runtime_error(
                    "DiscreteSparseBoundaryOperator::asMatrix(): "
                    "Epetra_CrsMatrix::ExtractMyRowView()) failed");
        if (transposed)
            for (int entry = 0; entry < entryCount; ++entry)
                mat(indices[entry], row) = values[entry];
        else
            for (int entry = 0; entry < entryCount; ++entry)
                mat(row, indices[entry]) = values[entry];
    }
    return mat;
}

template <typename ValueType>
unsigned int DiscreteSparseBoundaryOperator<ValueType>::rowCount() const
{
    return isTransposed() ? m_mat->NumGlobalCols() : m_mat->NumGlobalRows();
}

template <typename ValueType>
unsigned int DiscreteSparseBoundaryOperator<ValueType>::columnCount() const
{
    return isTransposed() ? m_mat->NumGlobalRows() : m_mat->NumGlobalCols();
}

template <typename ValueType>
void DiscreteSparseBoundaryOperator<ValueType>::addBlock(
        const std::vector<int>& rows,
        const std::vector<int>& cols,
        const ValueType alpha,
        arma::Mat<ValueType>& block) const
{
    // indices of entries of the untransposed (stored) matrix
    bool transposed = isTransposed();
    const std::vector<int>& untransposedRows = transposed ? cols : rows;
    const std::vector<int>& untransposedCols = transposed ? rows : cols;

    if (block.n_rows != rows.size() || block.n_cols != cols.size())
        throw std::invalid_argument(
                "DiscreteSparseBoundaryOperator::addBlock(): "
                "incorrect block size");

    int entryCount = 0;
    double* values = 0;
    int* indices = 0;

    for (size_t row = 0; row < untransposedRows.size(); ++row)
    {
        // Provision for future MPI support.
        if (m_mat->IndicesAreLocal())
        {
            int errorCode = m_mat->ExtractMyRowView(
                        untransposedRows[row], entryCount, values, indices);
            if (errorCode != 0)
                throw std::runtime_error(
                        "DiscreteSparseBoundaryOperator::addBlock(): "
                        "Epetra_CrsMatrix::ExtractMyRowView()) failed");
        }
        else
        {
            int errorCode = m_mat->ExtractGlobalRowView(
                        untransposedRows[row], entryCount, values, indices);
            if (errorCode != 0)
                throw std::runtime_error(
                        "DiscreteSparseBoundaryOperator::addBlock(): "
                        "Epetra_CrsMatrix::ExtractGlobalRowView()) failed");
        }

        for (size_t col = 0; col < untransposedCols.size(); ++col)
            for (int entry = 0; entry < entryCount; ++entry)
                if (indices[entry] == cols[col])
                    block(transposed ? col : row, transposed ? row : col) +=
                            alpha * static_cast<ValueType>(values[entry]);
    }
}

template <typename ValueType>
shared_ptr<const DiscreteSparseBoundaryOperator<ValueType> >
DiscreteSparseBoundaryOperator<ValueType>::castToSparse(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >&
        discreteOperator)
{
    shared_ptr<const DiscreteSparseBoundaryOperator<ValueType> > result =
        boost::dynamic_pointer_cast<const DiscreteSparseBoundaryOperator<ValueType> >(
                discreteOperator);
    if (result.get()==0 && (discreteOperator.get()!=0)) throw std::bad_cast();
    return result;

}


template <typename ValueType>
shared_ptr<const Epetra_CrsMatrix>
DiscreteSparseBoundaryOperator<ValueType>::epetraMatrix() const
{
    return m_mat;
}

template <typename ValueType>
TranspositionMode
DiscreteSparseBoundaryOperator<ValueType>::transpositionMode() const
{
    return m_trans;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteSparseBoundaryOperator<ValueType>::domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteSparseBoundaryOperator<ValueType>::range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool DiscreteSparseBoundaryOperator<ValueType>::opSupportedImpl(
        Thyra::EOpTransp M_trans) const
{
    return (M_trans == Thyra::NOTRANS || M_trans == Thyra::TRANS ||
            M_trans == Thyra::CONJ || M_trans == Thyra::CONJTRANS);
}

template <typename ValueType>
bool DiscreteSparseBoundaryOperator<ValueType>::isTransposed() const
{
    return m_trans & (TRANSPOSE | CONJUGATE_TRANSPOSE);
}

template <typename ValueType>
void DiscreteSparseBoundaryOperator<ValueType>::applyBuiltInImpl(
        const TranspositionMode trans,
        const arma::Col<ValueType>& x_in,
        arma::Col<ValueType>& y_inout,
        const ValueType alpha,
        const ValueType beta) const
{    
    TranspositionMode realTrans = trans;
    bool transposed = isTransposed();
    if (transposed)
        switch (trans) {
        case NO_TRANSPOSE: realTrans = TRANSPOSE; break;
        case TRANSPOSE: realTrans = NO_TRANSPOSE; break;
        case CONJUGATE: realTrans = CONJUGATE_TRANSPOSE; break;
        case CONJUGATE_TRANSPOSE: realTrans = CONJUGATE; break;
        // default: should not happen; anyway, don't change trans
        }

    reallyApplyBuiltInImpl(*m_mat, realTrans, x_in, y_inout, alpha, beta);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteSparseBoundaryOperator);

} // namespace Bempp

#endif // WITH_TRILINOS
