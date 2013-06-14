// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_LOCALGEOMETRY_HH
#define DUNE_FOAMGRID_LOCALGEOMETRY_HH

/**
 * \file
 * \brief The local geometry of th elements of FoamGrid
 */

namespace Dune{
    
    template<int mydim, int coorddim, class GridImp>
    class FoamGridLocalGeometry :
        public GenericGeometry::BasicGeometry<mydim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,mydim,coorddim> >
    {
        typedef GenericGeometry
        ::BasicGeometry<2, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,
                                                                  mydim,
                                                                  coorddim> > Base;
        
    public:
        /**
         * \brief Empty default constructor.
         */
        FoamGridLocalGeometry() 
        {}

        /**
         * \brief Construct geometry from coordinate vector.
         * \param type The type of the geometry.
         * \param coordinates The vector with the coordinates of the corners within
         * the father element.
         */
        FoamGridLocalGeometry(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,coorddim> >& coordinates) :
        Base(type, coordinates)
    {} 
    };
} // namespace Dune

#endif
