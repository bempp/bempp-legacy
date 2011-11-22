#ifndef DUNE_FOAMGRID_FACTORY_HH
#define DUNE_FOAMGRID_FACTORY_HH

/** \file
    \brief The specialization of the generic GridFactory for FoamGrid
    \author Oliver Sander
 */

#include <vector>
#include <map>

#include <dune/common/fvector.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

namespace Dune {

    /** \brief Specialization of the generic GridFactory for FoamGrid
        
    */
    template <int dimworld>
    class GridFactory<FoamGrid<dimworld> > 
        : public GridFactoryInterface<FoamGrid<dimworld> > {

        /** \brief Type used by the grid for coordinates */
        typedef typename FoamGrid<dimworld>::ctype ctype;

        typedef typename std::map<FieldVector<ctype,1>, unsigned int>::iterator VertexIterator;
        
        enum {dim = FoamGrid<dimworld>::dimension};

    public:

        /** \brief Default constructor */
        GridFactory()
            : factoryOwnsGrid_(true),
              vertexIndex_(0)
        {
            grid_ = new FoamGrid<dimworld>;

            grid_->entityImps_.resize(1);
        }

        /** \brief Constructor for a given grid object 

        If you already have your grid object constructed you can
        hand it over using this constructor.

        If you construct your factory class using this constructor
        the pointer handed over to you by the method createGrid() is
        the one you supplied here.
         */
        GridFactory(FoamGrid<dimworld>* grid)
            : factoryOwnsGrid_(false),
              vertexIndex_(0)
        {
            grid_ = grid;

            grid_->entityImps_.resize(1);
        }
        
        /** \brief Destructor */
        virtual ~GridFactory() {
            if (grid_ && factoryOwnsGrid_)
                delete grid_;
        }

        /** \brief Insert a vertex into the coarse grid */
        virtual void insertVertex(const FieldVector<ctype,dimworld>& pos) {
            Dune::get<0>(grid_->entityImps_[0]).push_back(FoamGridEntityImp<0,dimworld> (0,   // level
                                                                         pos,  // position
                                                                         grid_->freeIdCounter_[0]++));
            vertexArray_.push_back(&*Dune::get<0>(grid_->entityImps_[0]).rbegin());
        }

        /** \brief Insert an element into the coarse grid
            \param type The GeometryType of the new element
            \param vertices The vertices of the new element, using the DUNE numbering
        */
        virtual void insertElement(const GeometryType& type,
                                   const std::vector<unsigned int>& vertices) {

            assert(type.isTriangle());

            FoamGridEntityImp<2,dimworld> newElement(0,   // level
                                       grid_->freeIdCounter_[dim]++);  // id
            newElement.vertex_[0] = vertexArray_[vertices[0]];
            newElement.vertex_[1] = vertexArray_[vertices[1]];
            newElement.vertex_[2] = vertexArray_[vertices[2]];

            Dune::get<dim>(grid_->entityImps_[0]).push_back(newElement);

            
        }

        /** \brief Insert a boundary segment.

        This is only needed if you want to control the numbering of the boundary segments
        */
        virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices) {
            DUNE_THROW(Dune::NotImplemented, "insertBoundarySegment not implemented yet!");
        }

        /** \brief Insert a boundary segment (== a line) and the boundary segment geometry
         * 
            This influences the ordering of the boundary segments.
            Currently, the BoundarySegment object does not actually have any effect.
        */
        virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                                           const shared_ptr<BoundarySegment<2, 3> > boundarySegment)
        {
            insertBoundarySegment(vertices);
        }

        /** \brief Finalize grid creation and hand over the grid

        The receiver takes responsibility of the memory allocated for the grid
        */
        virtual FoamGrid<dimworld>* createGrid() {
            // Prevent a crash when this method is called twice in a row
            // You never know who may do this...
            if (grid_==NULL)
                return NULL;

            // ////////////////////////////////////////////////////
            //   Create the edges
            // ////////////////////////////////////////////////////

            // for convenience
            typedef FoamGridEntityImp<0,dimworld> FoamGridVertex;

            // For fast retrieval: a map from pairs of vertices to the edge that connects them
            std::map<std::pair<const FoamGridEntityImp<0,dimworld>*, const FoamGridEntityImp<0,dimworld>*>, FoamGridEntityImp<1,dimworld>*> edgeMap;

            typename std::list<FoamGridEntityImp<2,dimworld> >::iterator eIt    = Dune::get<2>(grid_->entityImps_[0]).begin();
            typename std::list<FoamGridEntityImp<2,dimworld> >::iterator eEndIt = Dune::get<2>(grid_->entityImps_[0]).end();

            for (; eIt!=eEndIt; ++eIt) {

                FoamGridEntityImp<2,dimworld>* element = &(*eIt);

                const Dune::GenericReferenceElement<double,dim>& refElement
                    = Dune::GenericReferenceElements<double, dim>::general(eIt->type());

                // Loop over all edges of this element
                for (size_t i=0; i<element->edges_.size(); ++i) {

                    // Get two vertices of the potential edge
                    const FoamGridVertex* v0 = element->vertex_[refElement.subEntity(i, 1, 0, 2)];
                    const FoamGridVertex* v1 = element->vertex_[refElement.subEntity(i, 1, 1, 2)];

                    FoamGridEntityImp<1,dimworld>* existingEdge = NULL;
                    typename std::map<std::pair<const FoamGridEntityImp<0,dimworld>*, const FoamGridEntityImp<0,dimworld>*>, FoamGridEntityImp<1,dimworld>*>::const_iterator e = edgeMap.find(std::make_pair(v0,v1));

                    if (e != edgeMap.end()) {
                        existingEdge = e->second;
                    } else {
                        e = edgeMap.find(std::make_pair(v1,v0));
                        if (e != edgeMap.end())
                            existingEdge = e->second;
                    }

                    if (existingEdge == NULL) {

                        // The current edge has not been inserted already.  We do that now

                        Dune::get<1>(grid_->entityImps_[0]).push_back(FoamGridEntityImp<1,dimworld>(v0,
                                                                                                    v1,
                                                                                                    0,   // level
                                                                                                    grid_->freeIdCounter_[1]++    // id
                                                                                                    ));

                        existingEdge = &*Dune::get<1>(grid_->entityImps_[0]).rbegin();
                                                                   
                        edgeMap.insert(std::make_pair(std::make_pair(v0,v1), existingEdge));
                        
                    }

                    // make element know about the edge
                    element->edges_[i] = existingEdge;

                    // make edge know about the element
                    existingEdge->elements_.push_back(element);

                }

            }

            // Create the index sets
            grid_->setIndices();
            

            // ////////////////////////////////////////////////
            //   Set the boundary ids
            //  \todo It must be possible to set them by hand
            // ////////////////////////////////////////////////

            unsigned int boundaryIdCounter = 0;

            for (typename std::list<FoamGridEntityImp<1,dimworld> >::iterator it = Dune::get<1>(grid_->entityImps_[0]).begin();
                 it != Dune::get<1>(grid_->entityImps_[0]).end();
                 ++it)
                it->boundaryId_ = boundaryIdCounter++;


            // ////////////////////////////////////////////////
            //   Hand over the new grid
            // ////////////////////////////////////////////////

            Dune::FoamGrid<dimworld>* tmp = grid_;
            grid_ = NULL;
            return tmp;
        }

    private:

        // Initialize the grid structure in UG
        void createBegin();

        // Pointer to the grid being built
        FoamGrid<dimworld>* grid_;

        // True if the factory allocated the grid itself, false if the
        // grid was handed over from the outside
        bool factoryOwnsGrid_;

        /** \brief Counter that creates the vertex indices */
        unsigned int vertexIndex_;

        std::vector<FoamGridEntityImp<0,dimworld>*> vertexArray_;

    };

}

#endif
