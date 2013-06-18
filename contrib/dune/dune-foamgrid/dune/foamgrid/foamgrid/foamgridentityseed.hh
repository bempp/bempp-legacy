#ifndef DUNE_FOAMGRID_ENTITY_SEED_HH
#define DUNE_FOAMGRID_ENTITY_SEED_HH

/**
 * \file
 * \brief The SubGridEntitySeed class
 */


namespace Dune {


/**
 * \brief The EntitySeed class provides the minmal information needed to restore an Entity using the grid.
 * \ingroup SubGrid
 *
 */
template<int codim, class GridImp>
class FoamGridEntitySeed
{
        template<int dimworld>
        friend class FoamGrid;

    protected:

        enum {dim = GridImp::dimension};
        enum {dimworld = GridImp::dimensionworld};
        enum {mydim = dim-codim};

        // Entity type of the hostgrid
        typedef FoamGridEntityImp<mydim, dimworld> EntityImplType;

    public:

        enum {codimension = codim};

        /**
         * \brief Create EntitySeed from hostgrid Entity
         *
         * We call hostEntity.seed() directly in the constructor
         * of SubGridEntitySeed to allow for return value optimization.
         *
         * If would use SubGridEntitySeed(hostEntity.seed())
         * we would have one copy even with optimization enabled.
         */
        FoamGridEntitySeed(const EntityImplType* impl) :
            entityImplPointer_(impl)
        {}

    protected:

        const EntityImplType* getImplementationPointer() const
        {
            return entityImplPointer_;
        }

    private:

        const FoamGridEntityImp<mydim, dimworld>* entityImplPointer_;
};

} // namespace Dune


#endif
