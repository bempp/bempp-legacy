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

class Geometry
{
public:
    virtual ~Geometry() {}
    virtual GeometryType type() const = 0;
    virtual bool affine() const = 0;
    virtual int cornerCount() const = 0;
    
    %apply arma::Mat<ctype>& ARGOUT_MAT { arma::Mat<ctype>& c };
    virtual void corners(arma::Mat<ctype>& c) const = 0;
    %clear arma::Mat<ctype>& c;

    %apply const arma::Mat<ctype>& IN_MAT { const arma::Mat<ctype>& local };
    %apply arma::Mat<ctype>& ARGOUT_MAT { arma::Mat<ctype>& global_ };
    virtual void local2global(const arma::Mat<ctype>& local,
                            arma::Mat<ctype>& global_) const = 0;
    %clear const arma::Mat<ctype>& local;
    %clear arma::Mat<ctype>& global_;

    %apply const arma::Mat<ctype>& IN_MAT { const arma::Mat<ctype>& global_ };
    %apply arma::Mat<ctype>& ARGOUT_MAT { arma::Mat<ctype>& local };
    virtual void global2local(const arma::Mat<ctype>& global_,
                            arma::Mat<ctype>& local) const = 0;
    %clear const arma::Mat<ctype>& global_;
    %clear arma::Mat<ctype>& local;

    %apply const arma::Mat<ctype>& IN_MAT { const arma::Mat<ctype>& local };
    %apply arma::Row<ctype>& ARGOUT_ROW { arma::Row<ctype>& int_element };
    virtual void integrationElement(const arma::Mat<ctype>& local,
                                    arma::Row<ctype>& int_element) const = 0;
    %clear const arma::Mat<ctype>& local;
    %clear arma::Row<ctype>& int_element;

    virtual ctype volume() const = 0;

    %apply arma::Col<ctype>& ARGOUT_COL { arma::Col<ctype>& c };
    virtual void center(arma::Col<ctype>& c) const = 0;
    %clear arma::Col<ctype>& c;

    %apply const arma::Mat<ctype>& IN_MAT { const arma::Mat<ctype>& local };
    %apply arma::Cube<ctype>& ARGOUT_CUBE { arma::Cube<ctype>& jacobian_t };
    virtual void jacobianTransposed(const arma::Mat<ctype>& local,
                                    arma::Cube<ctype>& jacobian_t) const = 0;
    %clear const arma::Mat<ctype>& local;
    %clear arma::Cube<ctype>& jacobian_t;

    %apply const arma::Mat<ctype>& IN_MAT { const arma::Mat<ctype>& local };
    %apply arma::Cube<ctype>& ARGOUT_CUBE { arma::Cube<ctype>& jacobian_inv_t };
    virtual void jacobianInverseTransposed(const arma::Mat<ctype>& local,
                                        arma::Cube<ctype>& jacobian_inv_t) const = 0;
    %clear const arma::Mat<ctype>& local;
    %clear arma::Cube<ctype>& jacobian_inv_t;
};

} // namespace Bempp
