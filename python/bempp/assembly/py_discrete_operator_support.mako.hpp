#ifndef BEMPP_PYTHON_DISCRETE_OPERATOR_SUPPORT_HPP
#define BEMPP_PYTHON_DISCRETE_OPERATOR_SUPPORT_HPP

#include "bempp/assembly/boundary_operator.hpp"
#include "bempp/utils/py_types.hpp"
#include "bempp/space/py_space_variants.hpp"
#include "bempp/assembly/discrete_sparse_boundary_operator.hpp"
#include "bempp/assembly/discrete_dense_boundary_operator.hpp"
#include "bempp/utils/py_types.hpp"
#include "fiber/scalar_traits.hpp"
#include <Python.h>
#include "numpy/arrayobject.h"
#include "Epetra_CrsMatrix.h"
#include <boost/variant.hpp>
#include <type_traits>


namespace Bempp {


template <typename ValueType>
PyObject* py_get_sparse_from_discrete_operator(const shared_ptr<const DiscreteBoundaryOperator< ValueType > >& op)
        {
                shared_ptr<const Epetra_CrsMatrix> epetraOperator = Bempp::DiscreteSparseBoundaryOperator< ValueType >::castToSparse(op)->epetraMatrix();

                PyObject* data;
                PyObject* colind;
                PyObject* row_ptr;
                PyObject* pModule;
                PyObject* pDict;
                PyObject* pFunc;
                PyObject* pArgsCsr;
                PyObject* pArgsShape;
                PyObject* pArgs;
                PyObject* pM;
                PyObject* pN;
                PyObject* pCsr;

                int M;
                int N;

                // Get the CSR matrix function

                pModule = PyImport_ImportModule("scipy.sparse");
                pDict = PyModule_GetDict(pModule);
                pFunc = PyDict_GetItemString(pDict, "csr_matrix");

                // Now get the parameter objects


                int* indexOffset;
                int* indices;
                double* values;

                npy_intp num_nonzeros = epetraOperator->NumGlobalNonzeros();
                npy_intp num_row_ptr = 1+epetraOperator->NumGlobalRows();

                M = epetraOperator->NumGlobalRows();
                N = epetraOperator->NumGlobalCols();

                epetraOperator->ExtractCrsDataPointers(indexOffset,indices,values);

                
                pM = PyInt_FromLong(M);
                pN = PyInt_FromLong(N);

                pArgsShape = PyTuple_New(2); 
                PyTuple_SetItem(pArgsShape, 0, pM);
                PyTuple_SetItem(pArgsShape, 1, pN);


                data = PyArray_SimpleNewFromData(1,&num_nonzeros,NPY_FLOAT64,values);
                colind = PyArray_SimpleNewFromData(1,&num_nonzeros,NPY_INT,indices);
                row_ptr = PyArray_SimpleNewFromData(1,&num_row_ptr,NPY_INT,indexOffset);
                
                
                pArgsCsr = PyTuple_New(3); 
                pArgs = PyTuple_New(2); 
                             
                PyTuple_SetItem(pArgsCsr, 0, data);
                PyTuple_SetItem(pArgsCsr, 1, colind);
                PyTuple_SetItem(pArgsCsr, 2, row_ptr);

                PyTuple_SetItem(pArgs,0,pArgsCsr);
                PyTuple_SetItem(pArgs,1,pArgsShape);
                
                // Now create the CSR matrix object
                
                pCsr = PyObject_CallObject(pFunc, pArgs);

                // Now clean up


                Py_DECREF(pModule);
                Py_DECREF(pArgs);

                return pCsr;

        }                


template <typename ValueType>
PyObject* py_array_from_dense_operator(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> op){

    return (*static_pointer_cast<const DiscreteDenseBoundaryOperator<ValueType>,
            const DiscreteBoundaryOperator<ValueType>>(op)).asNumpyObject();
}


}




#endif
