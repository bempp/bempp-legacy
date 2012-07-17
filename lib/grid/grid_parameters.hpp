#ifndef grid_parameters_hpp
#define grid_parameters_hpp

namespace Bempp {

/** \brief %Grid parameters.

  This structure is used to specify parameters of grid constructed by GridFactory.
  */
struct GridParameters {
    /** \brief %Grid topology */
    enum Topology {
        /** \brief one-dimensional grid embedded in a two-dimensional space */
        LINEAR,
        /** \brief two-dimensional grid composed of triangular elements,
            embedded in a three-dimensional space */
        TRIANGULAR,
        /** \brief two-dimensional grid composed of quadrilateral elements,
            embedded in a three-dimensional space */
        QUADRILATERAL,
        /** \brief two-dimensional grid composed (potentially) both of
            triangular and quadrilateral elements, embedded in a
            three-dimensional space */
        HYBRID_2D,
        /** \brief three-dimensional grid composed of tetrahedral elements,
            embedded in a three-dimensional space*/
        TETRAHEDRAL
    } topology;
};

}
#endif // grid_parameters_hpp
