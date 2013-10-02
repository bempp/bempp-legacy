%{
#include "assembly/discrete_boundary_operator.hpp"
#include <iostream>
#include <complex>
#include <armadillo>
%}

// TODO
// %include "discrete_boundary_operator_docstrings.i"

// For DiscreteBoundaryOperator, ignore the following warnings:
//   Nothing known about base class 'Thyra::LinearOpDefaultBase< ... >'
//   Nothing known about base class 'boost::enable_shared_from_this<
//     Bempp::DiscreteBoundaryOperator< ... > >'
//   Nothing known about 'Thyra::LinearOpDefaultBase< ... >::apply
// Reason: we don't need to create wrappers for these base classes,
// and providing their dummy definitions for SWIG
// leads to compilation errors for external modules using "%import bempp.swg"
// (because when %import is used, %{ ... %} blocks are ignored, so the
// wrapper code does not have the necessary includes).
%warnfilter(315,401) DiscreteBoundaryOperator< float >;
%warnfilter(315,401) DiscreteBoundaryOperator< double >;
%warnfilter(315,401) DiscreteBoundaryOperator< std::complex<float> >;
%warnfilter(315,401) DiscreteBoundaryOperator< std::complex<double> >;

%shared_ptr(Bempp::DiscreteBoundaryOperator<float>);
%shared_ptr(Bempp::DiscreteBoundaryOperator<double>);
%shared_ptr(Bempp::DiscreteBoundaryOperator<std::complex<float> >);
%shared_ptr(Bempp::DiscreteBoundaryOperator<std::complex<double> >);

#define shared_ptr boost::shared_ptr
namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(DiscreteBoundaryOperator);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_VALUE(DiscreteBoundaryOperator)

%extend DiscreteBoundaryOperator
{
    // this function is only for internal use
    %ignore addBlock;

    %apply const arma::Mat<float>& IN_MAT {
        const arma::Mat<float>& x_in
    };
    %apply const arma::Mat<double>& IN_MAT {
        const arma::Mat<double>& x_in
    };
    %apply const arma::Mat<std::complex<float> >& IN_MAT {
        const arma::Mat<std::complex<float> >& x_in
    };
    %apply const arma::Mat<std::complex<double> >& IN_MAT {
        const arma::Mat<std::complex<double> >& x_in
    };

    %apply arma::Mat<float>& INPLACE_MAT {
        arma::Mat<float>& y_inout
    };
    %apply arma::Mat<double>& INPLACE_MAT {
        arma::Mat<double>& y_inout
    };
    %apply arma::Mat<std::complex<float> >& INPLACE_MAT {
        arma::Mat<std::complex<float> >& y_inout
    };
    %apply arma::Mat<std::complex<double> >& INPLACE_MAT {
        arma::Mat<std::complex<double> >& y_inout
    };

    %apply arma::Mat<float>& ARGOUT_MAT {
        arma::Mat<float>& mat_out
    };
    %apply arma::Mat<double>& ARGOUT_MAT {
        arma::Mat<double>& mat_out
    };
    %apply arma::Mat<std::complex<float> >& ARGOUT_MAT {
        arma::Mat<std::complex<float> >& mat_out
    };
    %apply arma::Mat<std::complex<double> >& ARGOUT_MAT {
        arma::Mat<std::complex<double> >& mat_out
    };

    %apply const arma::Mat<float>& IN_MAT {
        const arma::Mat<float >& mat_in
    };
    %apply const arma::Mat<double>& IN_MAT {
        const arma::Mat<double >& mat_in
    };
    %apply const arma::Mat<std::complex<float> >& IN_MAT {
        const arma::Mat<std::complex<float> >& mat_in
    };
    %apply const arma::Mat<std::complex<double> >& IN_MAT {
        const arma::Mat<std::complex<double> >& mat_in
    };


    void asMatrix(arma::Mat<ValueType>& mat_out)
    {
        mat_out = $self->asMatrix();
    }

    %ignore asMatrix;

    static shared_ptr<const DiscreteBoundaryOperator<ValueType> >
    __addImpl(const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1, const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2){
        return op1+op2;
    }

    static shared_ptr<const DiscreteBoundaryOperator<ValueType> >
    __subImpl(const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1, const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2){
        return op1-op2;
    }


    static shared_ptr<const DiscreteBoundaryOperator<ValueType> >
    __scalarMultImpl(ValueType scalar, const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op){
        return scalar*op;
    }

    static shared_ptr<const DiscreteBoundaryOperator<ValueType> >
    __opCompositionImpl(const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
                        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2){
        return op1*op2;
    }

    static void
    __matrixMultImpl(const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
                     const arma::Mat<ValueType>& mat_in,
                     arma::Mat<ValueType>& mat_out){

        const int opRows = op->rowCount();
        const int opCols = op->columnCount();
        const int matCols = mat_in.n_cols;
        const int matRows = mat_in.n_rows;
        mat_out.zeros(opRows,matCols);
        if (opCols!=matRows)
            throw std::invalid_argument("__matrixMultImpl(): Wrong dimensions.");
        for (int i=0;i<matCols;i++){
            arma::Col< ValueType > tmp1 = mat_in.unsafe_col(i);
            arma::Col< ValueType > tmp2 = mat_out.unsafe_col(i);
            op->apply(Bempp::NO_TRANSPOSE,tmp1,tmp2,1.0,0.0);
        }
    }

    static void
    __matrixHMultImpl(const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
                     const arma::Mat<ValueType>& mat_in,
                     arma::Mat<ValueType>& mat_out){

        const int opRows = op->rowCount();
        const int opCols = op->columnCount();
        const int matCols = mat_in.n_cols;
        const int matRows = mat_in.n_rows;
        mat_out.zeros(opCols,matCols);
        if (opRows!=matRows)
            throw std::invalid_argument("__matrixHMultImpl(): Wrong dimensions.");
        for (int i=0;i<matCols;i++){
            arma::Col< ValueType > tmp1 = mat_in.unsafe_col(i);
            arma::Col< ValueType > tmp2 = mat_out.unsafe_col(i);
            op->apply(Bempp::CONJUGATE_TRANSPOSE,tmp1,tmp2,1.0,0.0);
        }
    }

    %feature("compactdefaultargs") asDiscreteAcaBoundaryOperator;

    %pythoncode {
        def __add__(self,other):
            return self.__addImpl(self,other)

        def __sub__(self,other):
            return self.__subImpl(self,other)

        def __mul__(self,other):
            from numbers import Number
            import numpy as np
            if isinstance(other,(Number,np.number)):
                return self.__scalarMultImpl(other,self)
            elif isinstance(other,np.ndarray):
                ndim = other.ndim
                if ndim==0:
                    return self.__scalarMultImpl(other,self)
                elif ndim==1:
                    return self.__matrixMultImpl(self,other[:,np.newaxis])[:,0]
                elif ndim==2:
                    return self.__matrixMultImpl(self,other)
                else:
                    raise ValueError("Discrete boundary operators do not support "
                                     "multiplication by arrays with more than 2 "
                                     "dimensions.")
            elif isinstance(other,type(self)):
                return self.__opCompositionImpl(self,other)
            else:
                raise ValueError("Discrete boundary operators do not support "
                                 "multiplication with this type.")

        def __div__(self,other):
            if other==0: raise ValueError("Division by zero not allowed.")
            return self.__mul__(1./other)

        def __rmul__(self,other):
            from numbers import Number
            import numpy as np
            if isinstance(other,(Number,np.number)):
                return self.__scalarMultImpl(other,self)
            elif isinstance(other,np.ndarray):
                if len(other.shape)==0:
                    return self.__scalarMultImpl(other,self)
                else:
                    raise ValueError("Discrete boundary operators do not support "
                                    "multiplication with this type.")
            else:
                raise ValueError("Discrete boundary operators do not support "
                                 "multiplication with this type.")

        def matvec(self, other):
            """
            Multiply this operator with 'other' and return the result.

            *Parameters:*
               - other (1D or 2D ndarray)
                    Vector or matrix that should be multiplied with this operator.
            """
            sn = len(other.shape)
            if sn == 1:
               data = other.reshape(other.shape[0],1)
            else:
               data = other
            res = self.__matrixMultImpl(self,data)
            if sn == 1:
               return res.reshape(self.rowCount())
            else:
               return res

        def rmatvec(self, other):
            """
            Multiply the conjugate transpose of this operator with 'other'
            and return the result.

            *Parameters:*
               - other (1D or 2D ndarray)
                    Vector or matrix that should be multiplied with the
                    conjugate transpose of this operator.
            """
            sn = len(other.shape)
            if sn == 1:
               data = other.reshape(other.shape[0],1)
            else:
               data = other
            res = self.__matrixHMultImpl(self,data)
            if sn == 1:
               return res.reshape(self.columnCount())
            else:
               return res

        def matmat(self, other):
            """
            Multiply this operator with the matrix 'other' and return the result.

            *Parameters:*
               - other (2D ndarray)
                    Matrix that should be multiplied with this operator.
            """
            return self.__matrixMultImpl(self,other)

        def asPyTrilinosOperator(self,label="",pid=0):
            """
            Return a representation as Trilinos Epetra Operator.

            *Parameters:*
               - label (str)
                    The operator's label.
               - pid (int)
                    Id of the process in which the operator lives (not yet supported).
            """
            import trilinos

            return trilinos.PyTrilinosOperator(self,label=label,pid=pid)

        def asPyTrilinosInverseOperator(self,label="",pid=0):
            """
            Return a representation as Trilinos Epetra Inverse Operator.

            *Parameters:*
               - label (str)
                    The operator's label.
               - pid (int)
                    Id of the process in which the operator lives (not yet supported).
            """
            import trilinos

            return trilinos.PyTrilinosInverseOperator(self,label=label,pid=pid)


        @property
        def shape(self):
            return (self.rowCount(),self.columnCount())

        @property
        def dtype(self):
            import numpy
            return numpy.dtype(self.valueType())
    }

} // %extend DiscreteBoundaryOperator

} // namespace Bempp

%include "assembly/discrete_boundary_operator.hpp"

#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(DiscreteBoundaryOperator);

%clear const arma::Mat<float>& x_in;
%clear const arma::Mat<double>& x_in;
%clear const arma::Mat<std::complex<float> >& x_in;
%clear const arma::Mat<std::complex<double> >& x_in;

%clear arma::Mat<float>& y_inout;
%clear arma::Mat<double>& y_inout;
%clear arma::Mat<std::complex<float> >& y_inout;
%clear arma::Mat<std::complex<double> >& y_inout;

%clear arma::Mat<float>& mat_out;
%clear arma::Mat<double>& mat_out;
%clear arma::Mat<std::complex<float> >& mat_out;
%clear arma::Mat<std::complex<double> >& mat_out;

%clear const arma::Mat<float>& mat_in;
%clear const arma::Mat<double>& mat_in;
%clear const arma::Mat<std::complex<float> >& mat_in;
%clear const arma::Mat<std::complex<double> >& mat_in;

} // namespace Bempp

