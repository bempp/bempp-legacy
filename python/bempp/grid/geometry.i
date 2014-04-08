%{
#include "grid/geometry.hpp"
%}

%include "geometry_docstrings.i"

/*
    This is fairly elegant but doesn't work yet because there are 
    both input and output arguments names "local" and "global",
    so SWIG confuses their typemaps.

    TODO: rename all output arguments to "outSomething".
*/
/*
namespace Bempp {
    %apply const arma::Mat<ctype>& IN_MAT {
        const arma::Mat<ctype>& local, 
        const arma::Mat<ctype>& global 
    };

    %apply arma::Col<ctype>& ARGOUT_COL {
        arma::Col<ctype>& c
    };

    %apply arma::Row<ctype>& ARGOUT_ROW {
        arma::Row<ctype>& int_element
    };

    %apply arma::Mat<ctype>& ARGOUT_MAT {
        arma::Mat<ctype>& c, 
        arma::Mat<ctype>& global, 
        arma::Mat<ctype>& local
    };

    %apply arma::Cube<ctype>& ARGOUT_CUBE {
        arma::Cube<ctype>& jacobian_t,
        arma::Cube<ctype>& jacobian_inv_t
    };
}

%include "grid/geometry.hpp"

namespace Bempp {
    %clear const arma::Mat<ctype>& local;
    %clear const arma::Mat<ctype>& global;

    %clear arma::Col<ctype>& c;

    %clear arma::Row<ctype>& int_element;

    %clear arma::Mat<ctype>& c;
    %clear arma::Mat<ctype>& global; 
    %clear arma::Mat<ctype>& local;

    %clear arma::Cube<ctype>& jacobian_t;
    %clear arma::Cube<ctype>& jacobian_inv_t;
}
*/

/* We resort to copying the class declaration... sigh. */

namespace Bempp
{

%extend Geometry 
{
    // Reference to the parent entity, stored to prevent it from
    // being garbage-collected while this geometry is alive
    %pythoncode %{
        def parentEntity(self): return self._parentEntity
        parentEntity = property(parentEntity, None, None, "Parent entity")
    %}
}

%rename(corners) Geometry::getCorners;
%rename(normals) Geometry::getNormals;
%rename(integrationElements) Geometry::getIntegrationElements;
%rename(center) Geometry::getCenter;
%rename(jacobiansTransposed) Geometry::getJacobiansTransposed;
%rename(jacobianInversesTransposed) Geometry::getJacobianInversesTransposed;

class Geometry
{
public:
    virtual ~Geometry() {}
    virtual GeometryType type() const = 0;
    virtual bool affine() const = 0;
    virtual int cornerCount() const = 0;
    
    %apply arma::Mat<ctype>& ARGOUT_MAT { arma::Mat<ctype>& c };
    virtual void getCorners(arma::Mat<ctype>& c) const = 0;
    %clear arma::Mat<ctype>& c;

    %apply const arma::Mat<ctype>& IN_MAT { const arma::Mat<ctype>& local };
    %apply const arma::Mat<ctype>& IN_MAT { const arma::Mat<ctype>& global_ };
    %apply arma::Mat<ctype>& ARGOUT_MAT { arma::Mat<ctype>& outLocal };
    %apply arma::Mat<ctype>& ARGOUT_MAT { arma::Mat<ctype>& outGlobal };
    %apply arma::Mat<ctype>& ARGOUT_MAT { arma::Mat<ctype>& outNormal };
#   ifdef NDEBUG
      %typemap(check) (const arma::Mat<ctype>& local,
                              arma::Mat<ctype>& outGlobal) {
          // $self does not seem to work here. So using arg1 explicitly...
          if((int)$1->n_rows != arg1->dim()) {
              PyErr_SetString(PyExc_RuntimeError, "Incorrect number of dimensions");
              SWIG_fail;
          }
      }
#   endif

    virtual void local2global(const arma::Mat<ctype>& local,
                            arma::Mat<ctype>& outGlobal) const = 0;
    virtual void global2local(const arma::Mat<ctype>& global_,
                            arma::Mat<ctype>& outLocal) const = 0;
    virtual void getNormals(const arma::Mat<ctype>& local,
                            arma::Mat<ctype>& outNormal) const = 0;

    %apply arma::Row<ctype>& ARGOUT_ROW { arma::Row<ctype>& outIntElement };
    virtual void getIntegrationElements(const arma::Mat<ctype>& local,
                                    arma::Row<ctype>& outIntElement) const = 0;
    %clear arma::Row<ctype>& outIntElement;

    virtual ctype volume() const = 0;

    %apply arma::Col<ctype>& ARGOUT_COL { arma::Col<ctype>& c };
    virtual void getCenter(arma::Col<ctype>& c) const = 0;
    %clear arma::Col<ctype>& c;

    %apply arma::Cube<ctype>& ARGOUT_CUBE { arma::Cube<ctype>& outJacobianT };
    virtual void getJacobiansTransposed(const arma::Mat<ctype>& local,
                                    arma::Cube<ctype>& outJacobianT) const = 0;
    %clear arma::Cube<ctype>& jacobian_t;

    %apply const arma::Mat<ctype>& IN_MAT { const arma::Mat<ctype>& local };
    %apply arma::Cube<ctype>& ARGOUT_CUBE { arma::Cube<ctype>& outJacobianInvT };
    virtual void getJacobianInversesTransposed(const arma::Mat<ctype>& local,
                                        arma::Cube<ctype>& outJacobianInvT) const = 0;
    %clear arma::Cube<ctype>& outJacobianInvT;

    %clear const arma::Mat<ctype>& local;
    %clear const arma::Mat<ctype>& global_;
    %clear arma::Mat<ctype>& outLocal;
    %clear arma::Mat<ctype>& outGlobal;
    %clear arma::Mat<ctype>& outNormal;
};

} // namespace Bempp
