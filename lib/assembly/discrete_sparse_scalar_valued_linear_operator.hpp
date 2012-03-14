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

#ifndef bempp_discrete_sparse_scalar_valued_linear_operator_hpp
#define bempp_discrete_sparse_scalar_valued_linear_operator_hpp

#include "discrete_scalar_valued_linear_operator.hpp"

#include <iostream>
#include <stdexcept>

#ifdef WITH_TRILINOS
#include <Epetra_FECrsMatrix.h>
#include <Epetra_SerialComm.h>
#endif

namespace Bempp {

template <typename ValueType>
class DiscreteSparseScalarValuedLinearOperator :
        public DiscreteScalarValuedLinearOperator<ValueType>
{
#ifdef WITH_TRILINOS
public:
    DiscreteSparseScalarValuedLinearOperator(std::auto_ptr<Epetra_FECrsMatrix> mat) :
        m_mat(mat) {}
#else
    // This class cannot be used without Trilinos
private:
    DiscreteSparseScalarValuedLinearOperator();
#endif

public:
    virtual void multiplyAddVector(ValueType multiplier,
                                   const arma::Col<ValueType>& argument,
                                   arma::Col<ValueType>& result)
    {
        throw std::runtime_error("DiscreteSparseScalarValuedLinearOperator::"
                                 "multiplyAddVector(): not implemented yet");
    }

    virtual void dump() const
    {
#ifdef WITH_TRILINOS
        std::cout << *m_mat << std::endl;
#endif
    }

    virtual arma::Mat<ValueType> asMatrix() const
    {
#ifdef WITH_TRILINOS
        if (m_mat->Comm().NumProc() != 1)
            throw std::runtime_error(
                    "DiscreteSparseScalarValuedLinearOperator::asMatrix(): "
                    "conversion of distributed matrices to local matrices is unsupported");

        const int rowCount = m_mat->NumGlobalRows();
        const int colCount = m_mat->NumGlobalCols();
        arma::Mat<ValueType> mat(rowCount, colCount);
        mat.fill(0.);
        for (int row = 0; row < rowCount; ++row)
        {
            int entryCount = 0;
            double* values = 0;
            int* indices = 0;
            int errorCode = m_mat->ExtractMyRowView(row, entryCount,
                                                    values, indices);
            if (errorCode != 0)
                throw std::runtime_error(
                        "DiscreteSparseScalarValuedLinearOperator::asMatrix(): "
                        "Epetra_CrsMatrix::ExtractMyRowView()) failed");
            for (int entry = 0; entry < entryCount; ++entry)
                mat(row, indices[entry]) = values[entry];
        }
        return mat;
#endif
    }

private:
#ifdef WITH_TRILINOS
    std::auto_ptr<Epetra_FECrsMatrix> m_mat;
#endif
};

} // namespace Bempp

#endif
