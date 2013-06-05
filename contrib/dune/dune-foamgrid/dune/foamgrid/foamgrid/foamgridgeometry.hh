#ifndef DUNE_FOAMGRID_GEOMETRY_HH
#define DUNE_FOAMGRID_GEOMETRY_HH

/** \file
* \brief The FoamGridGeometry class and its specializations
*/

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/geometry/genericgeometry/geometry.hh>


namespace Dune {

template<int mydim, int coorddim, class GridImp> class FoamGridGeometry;

namespace FacadeOptions
{
    template< int mydim, int cdim, class GridImp>
    struct StoreGeometryReference<mydim, cdim, GridImp, FoamGridGeometry>
    {
        static const bool v = false;
    };
}

template<int mydim, int coorddim, class GridImp>
class FoamGridGeometry :
        public GenericGeometry::BasicGeometry<mydim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,mydim,coorddim> >
{

    typedef typename GenericGeometry::BasicGeometry<mydim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,mydim,coorddim> > Base;

    public:

    /**
     * \brief This is DefaultConstructor
     */
    FoamGridGeometry() {}

    /**
     * \brief Construct geometry from coordinate vector
     */
    FoamGridGeometry(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,coorddim> >& coordinates) :
        Base(type, coordinates)
    {}

};


}  // namespace Dune

#endif
