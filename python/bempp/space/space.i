%{
#include "space/space.hpp"
#include <complex>
%}

%shared_ptr(Bempp::Space<float>);
%shared_ptr(Bempp::Space<double>);
%shared_ptr(Bempp::Space<std::complex<float> >);
%shared_ptr(Bempp::Space<std::complex<double> >);

// TODO
// %include "space_docstrings.i"

namespace Bempp
{

template<typename BasisFunctionType> class Space;

%extend Space
{

// this function is only for internal use
%ignore basis;
%ignore shapeset;

// this function is only for internal use
%ignore shapeFunctionValue;
%ignore basisFunctionValue;

// to be wrapped later...
%ignore setElementVariant;

// to be wrapped later...
%ignore elementVariant;

%ignore operator=;

%ignore getGlobalDofs;

%ignore global2localDofs;
%ignore flatLocal2localDofs;

// these functions are only for internal use
%ignore getGlobalDofPositions;
%ignore getFlatLocalDofPositions;
%ignore getGlobalDofBoundingBoxes;
%ignore getFlatLocalDofBoundingBoxes;
%ignore getGlobalDofNormals;
%ignore getFlatLocalDofNormals;

%ignore dumpClusterIds;
%ignore dumpClusterIdsEx;
}

%define BEMPP_EXTEND_SPACE(BASIS, PYBASIS)
    %extend Space< BASIS >
    {
        %pythonprepend assignDofs %{
            print ("HINT: It is not necessary to call Space.assignDofs() any more.\n"
                   "DOFs are now automatically assigned as soon as a space object is "
                   "constructed.")
        %}

        %apply arma::Mat<float>& ARGOUT_MAT {
            arma::Mat<float>& positions,
            arma::Mat<float>& normals
        };
        %apply arma::Mat<double>& ARGOUT_MAT {
            arma::Mat<double>& positions,
            arma::Mat<double>& normals
        };

        void globalDofPositions(arma::Mat<Bempp::Space< BASIS >::CoordinateType>& positions)
        {
            typedef Bempp::Space< BASIS >::CoordinateType CoordinateType;
            std::vector<Bempp::Point3D<CoordinateType> > vecPositions;
            $self->getGlobalDofPositions(vecPositions);
            const int dimWorld = $self->grid()->dimWorld();
            positions.set_size(dimWorld, vecPositions.size());
            for (size_t p = 0; p < vecPositions.size(); ++p) {
                positions(0, p) = vecPositions[p].x;
                if (dimWorld >= 2)
                    positions(1, p) = vecPositions[p].y;
                if (dimWorld >= 3)
                    positions(2, p) = vecPositions[p].z;
            }
        }

        void globalDofNormals(arma::Mat<Bempp::Space< BASIS >::CoordinateType>& normals)
        {
            typedef Bempp::Space< BASIS >::CoordinateType CoordinateType;
            std::vector<Bempp::Point3D<CoordinateType> > vecNormals;
            $self->getGlobalDofNormals(vecNormals);
            const int dimWorld = $self->grid()->dimWorld();
            normals.set_size(dimWorld, vecNormals.size());
            for (size_t p = 0; p < vecNormals.size(); ++p) {
                normals(0, p) = vecNormals[p].x;
                if (dimWorld >= 2)
                    normals(1, p) = vecNormals[p].y;
                if (dimWorld >= 3)
                    normals(2, p) = vecNormals[p].z;
            }
        }

    }
%enddef
BEMPP_ITERATE_OVER_BASIS_TYPES(BEMPP_EXTEND_SPACE);

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS(Space);

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "space/space.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(Space);

%clear const arma::Mat<float>& positions;
%clear const arma::Mat<float>& normals;
%clear const arma::Mat<double>& normals;
%clear const arma::Mat<double>& positions;
}

