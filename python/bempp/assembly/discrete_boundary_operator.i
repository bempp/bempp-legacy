%{
#include "assembly/discrete_boundary_operator.hpp"
#include <iostream>
#include <complex>
#include <armadillo>
%}

// TODO
// %include "discrete_boundary_operator_docstrings.i"

%shared_ptr(boost::enable_shared_from_this<Bempp::DiscreteBoundaryOperator<float> >)
%shared_ptr(boost::enable_shared_from_this<Bempp::DiscreteBoundaryOperator<double> >)
%shared_ptr(boost::enable_shared_from_this<Bempp::DiscreteBoundaryOperator<std::complex<float> > >)
%shared_ptr(boost::enable_shared_from_this<Bempp::DiscreteBoundaryOperator<std::complex<double> > >)

%shared_ptr(Thyra::LinearOpDefaultBase<float>);
%shared_ptr(Thyra::LinearOpDefaultBase<double>);
%shared_ptr(Thyra::LinearOpDefaultBase<std::complex<float> >);
%shared_ptr(Thyra::LinearOpDefaultBase<std::complex<double> >);

%shared_ptr(Bempp::DiscreteBoundaryOperator<float>);
%shared_ptr(Bempp::DiscreteBoundaryOperator<double>);
%shared_ptr(Bempp::DiscreteBoundaryOperator<std::complex<float> >);
%shared_ptr(Bempp::DiscreteBoundaryOperator<std::complex<double> >);


#define shared_ptr boost::shared_ptr
namespace Thyra
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(LinearOpDefaultBase);

template <typename ValueType>
class LinearOpDefaultBase
{
public:
    virtual ~LinearOpDefaultBase() = 0; // prevent instantiation
};

} // namespace Thyra

namespace boost
{

template <typename T> class enable_shared_from_this;

template <typename T>
class enable_shared_from_this
{
public:
    virtual ~enable_shared_from_this() = 0;
};

%template(enable_shared_from_this_discrete_boundary_operator_float32) enable_shared_from_this<Bempp::DiscreteBoundaryOperator<float> >;
%template(enable_shared_from_this_discrete_boundary_operator_float64) enable_shared_from_this<Bempp::DiscreteBoundaryOperator<double> >;
%template(enable_shared_from_this_discrete_boundary_operator_complex64) enable_shared_from_this<Bempp::DiscreteBoundaryOperator<std::complex<float> > >;
%template(enable_shared_from_this_discrete_boundary_operator_complex128) enable_shared_from_this<Bempp::DiscreteBoundaryOperator<std::complex<double> > >;

} // namespace boost



namespace Bempp
{

/* DECLARE_TEMPLATE_VALUE_METHOD_AUTO_DOCSTRING(DiscreteBoundaryOperator, apply); */

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(DiscreteBoundaryOperator);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_VALUE(DiscreteBoundaryOperator)

// Nothing known about 'Thyra::LinearOpDefaultBase< ... >::apply
%warnfilter(315) DiscreteBoundaryOperator;

%extend DiscreteBoundaryOperator
{
    // this function is only for internal use
    %ignore addBlock;

    %apply const arma::Col<float>& IN_COL {
        const arma::Col<float>& x_in
    };
    %apply const arma::Col<double>& IN_COL {
        const arma::Col<double>& x_in
    };
    %apply const arma::Col<std::complex<float> >& IN_COL {
        const arma::Col<std::complex<float> >& x_in
    };
    %apply const arma::Col<std::complex<double> >& IN_COL {
        const arma::Col<std::complex<double> >& x_in
    };

    %apply arma::Col<float>& INPLACE_COL {
        arma::Col<float>& y_inout
    };
    %apply arma::Col<double>& INPLACE_COL {
        arma::Col<double>& y_inout
    };
    %apply arma::Col<std::complex<float> >& INPLACE_COL {
        arma::Col<std::complex<float> >& y_inout
    };
    %apply arma::Col<std::complex<double> >& INPLACE_COL {
        arma::Col<std::complex<double> >& y_inout
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

    static void
    __matrixMultImpl(const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
                     const arma::Mat<ValueType>& mat_in,
                     arma::Mat<ValueType>& mat_out){

        const int opRows = op->rowCount();
        const int opCols = op->columnCount();
        const int matCols = mat_in.n_cols;
        const int matRows = mat_in.n_rows;
        mat_out.zeros(opRows,matCols);
        if (opCols!=matRows) throw std::runtime_error("__matrixMultImpl(op,mat_in,mat_out: Wrong dimensions.");
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
        if (opRows!=matRows) throw std::runtime_error("__matrixHMultImpl(op,mat_in,mat_out: Wrong dimensions.");
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
                if len(other.shape)==0:
                    return self.__scalarMultImpl(other,self)
                else:
                    return self.__matrixMultImpl(self,other)
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

        def matvec(self,other):
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
			
        def rmatvec(self,other):
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

        def matmat(self,other):
            return self.__matrixMultImpl(self,other)
	            

        @property
        def shape(self):
            return (self.rowCount(),self.columnCount())
    }



}

} // namespace Bempp

%include "assembly/discrete_boundary_operator.hpp"

#undef shared_ptr

namespace Thyra
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(LinearOpDefaultBase);
}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(DiscreteBoundaryOperator);

%clear const arma::Col<float>& x_in;
%clear const arma::Col<double>& x_in;
%clear const arma::Col<std::complex<float> >& x_in;
%clear const arma::Col<std::complex<double> >& x_in;

%clear arma::Col<float>& y_inout;
%clear arma::Col<double>& y_inout;
%clear arma::Col<std::complex<float> >& y_inout;
%clear arma::Col<std::complex<double> >& y_inout;

%clear arma::Mat<float>& mat_out;
%clear arma::Mat<double>& mat_out;
%clear arma::Mat<std::complex<float> >& mat_out;
%clear arma::Mat<std::complex<double> >& mat_out;

%clear const arma::Mat<float>& mat_in;
%clear const arma::Mat<double>& mat_in;
%clear const arma::Mat<std::complex<float> >& mat_in;
%clear const arma::Mat<std::complex<double> >& mat_in;


} // namespace Bempp

