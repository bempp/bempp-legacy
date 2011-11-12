// Copyright (C) 2011 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef bempp_lib_grid_3d_entity_decl_hpp
#define bempp_lib_grid_3d_entity_decl_hpp

#include "geometry_decl.hpp"

#include <dune/common/static_assert.hh>
#include <dune/grid/common/geometry.hh>

namespace Bempp {

// Forward declarations
class Geometry;
template<int codim> class EntityPointer;
template<int codim> class EntityIterator;

/** \brief Abstract wrapper of an entity.

 \param codim entity's codimension
 */
template<int codim>
class Entity {
public:
	/** Destructor */
	virtual ~Entity() {
	}

	/** Entity codimension */
	enum {
		codimension = codim
	};

	/** Entity level */
	virtual int level() const = 0;

	/** \brief Reference to the geometry of this entity.

	 This object gives, among other things, the map from a reference element to world coordinates.

	 \note Be careful when storing such references. If the state
	 of any object is changed, e.g. an iterator is advanced, there
	 is no guarantee that the reference remains valid.
	 */
	virtual const Geometry& geometry() const = 0;

	/** \brief Return the name of the reference element. The type can
	 be used to access the Dune::GenericReferenceElement.
	 */
	virtual Dune::GeometryType type() const = 0;
};

/** Abstract wrapper of an entity of codimension 0. */
template<>
class Entity<0> {
public:
	/** Destructor */
	virtual ~Entity() {
	}

	/** Entity codimension */
	enum {
		codimension = 0
	};

	/** @name Methods shared by entities of all codimensions @{ */

	/** Entity level */
	virtual int level() const = 0;

	/** \brief Reference to the geometry of this entity.

	 This object gives, among other things, the map from a reference element to world coordinates.

	 \note Be careful when storing such references. If the state
	 of any object is changed, e.g. an iterator is advanced, there
	 is no guarantee that the reference remains valid.
	 */
	virtual const Geometry& geometry() const = 0;

	/** \brief Return the name of the reference element. The type can
	 be used to access the Dune::GenericReferenceElement.
	 */
	virtual Dune::GeometryType type() const = 0;

	/** @} @name Extended interface of entities of codimension 0 @{ */

	/** \brief Number of vertices. This method is in principle
	 redundant because this information can be obtained via the
	 reference element of the geometry. It is there for efficiency
	 reasons and to make the interface self-contained.
	 */
	virtual int vertexCount() const = 0;

	/** \brief Number of edges. This method is in principle
	 redundant because this information can be obtained via the
	 reference element of the geometry. It is there for efficiency
	 reasons and to make the interface self-contained. */
	virtual int edgeCount() const = 0;

	/** \brief Iterator over subedges.

	 The user is responsible for deleting the returned iterator */
	virtual EntityIterator<1>* subEdgeIterator() const = 0;

	/** \brief Iterator over subvertices.

	 The user is responsible for deleting the returned iterator */
	virtual EntityIterator<2>* subVertexIterator() const = 0;

	// To be implemented later.
	//  /** \brief Access to intersections with neighboring leaf elements.
	//     A neighbor is an entity of codimension 0
	//     which has an intersection of codimension 1 in common with this entity.
	//     Access to those neighbors is provided using the IntersectionIterator.
	//     This method returns an iterator refering to the first neighbor.

	//     \note If the partitionType of the Entity is GhostEntity,
	//     this method might give you only one neighbor, which is the
	//     interior Entity the GhostEntity is connected to.
	//  */
	//  virtual SurfaceGridIntersectionIterator ileafbegin() const = 0;

	//  /** \brief Reference to an IntersectionIterator one past the last intersection

	//     \note If the partitionType of the Entity is GhostEntity,
	//     this method might give you only one neighbor, which is the
	//     interior Entity the GhostEntity is connected to.
	//  */
	//  virtual SurfaceGridIntersectionIterator ileafend() const = 0;

	//  /** \brief Intra-level access to intersections with neighboring elements.
	//     A neighbor is an entity of codimension 0
	//     which has an intersection of codimension 1 in common with this entity.
	//     Access to those neighbors is provided using the IntersectionIterator.
	//     This method returns an iterator refering to the first neighbor.

	//     \note If the partitionType of the Entity is GhostEntity,
	//     this method might give you only one neighbor, which is the
	//     interior Entity the GhostEntity is connected to.
	//  */
	//  virtual SurfaceGridIntersectionIterator ilevelbegin() const = 0;

	//  /** \brief Reference to an IntersectionIterator one past the last intersection

	//     \note If the partitionType of the Entity is GhostEntity,
	//     this method might give you only one neighbor, which is the
	//     interior Entity the GhostEntity is connected to.
	//  */
	//  virtual SurfaceGridIntersectionIterator ilevelend() const = 0;

	/** \brief Inter-level access to father entity on the next-coarser grid.
	 The given entity resulted directly from a subdivision of its father
	 entity. For the macro elements dereferencing the EntityPointer is undefined.

	 The caller is responsible for deleting the returned entity pointer.

	 \note If the partitionType of the Entity is GhostEntity,
	 it is not guaranteed that this method is working
	 or implemented in general.
	 For some grids it might be available, though.
	 */
	virtual EntityPointer<0>* father() const = 0;

	/** \brief Return true if entity has a father entity which can be accessed
	 using the father() method.
	 */
	virtual bool hasFather() const = 0;

	/** Returns true if the entity is contained in the leaf grid */
	virtual bool isLeaf() const = 0;

	/** \brief Returns true if element is of regular type in red/green type refinement.
	 In bisection or hanging node refinement this is always true.
	 */
	virtual bool isRegular() const = 0;

	// To be implemented later.
	//  /** \brief Provides information how this element has been subdivided from
	//     its father element.
	//     The returned LocalGeometry is a model of Dune::Geometry<dimension,dimension,...>
	//     mapping from the reference element of the given element to the reference
	//     element of the father element.
	//     This is sufficient to interpolate all degrees of freedom in the
	//     conforming case. Nonconforming may require access to neighbors of father and
	//     computations with local coordinates.
	//     On the fly case is somewhat inefficient since degrees of freedom
	//     may be visited several times.
	//     If we store interpolation matrices, this is tolerable. We assume that on-the-fly
	//     implementation of interpolation is only done for simple discretizations.

	//     \note If the partitionType of the Entity is GhostEntity,
	//     it is not guaranteed that this method is working
	//     or implemented in general.
	//     For some grids it might be available, though.
	//  */
	//  virtual const SurfaceGridLocalGeometry& geometryInFather() const = 0;

	/** \brief Inter-level access to elements that resulted from (recursive)
	 subdivision of this element.

	 \param[in] maxlevel Iterator does not stop at elements with level greater than maxlevel.
	 \return Iterator to the first son (level is not greater than maxlevel)

	 The caller is responsible for deleting the returned iterator.

	 \note If the partitionType of the Entity is GhostEntity,
	 it is not guaranteed that this method is working
	 or implemented in general.
	 For some grids it might be available, though.
	 */
	virtual EntityIterator<0>* sonIterator(int maxlevel) const = 0;

	/** \brief Returns true if the entity has been created during the last call to adapt()
	 */
	virtual bool isNew() const = 0;

	/** \brief Returns true if the entity might disappear during the next call to adapt().
	 * If the method returns false, the entity is guaranteed to still be present after
	 * adaptation.
	 */
	virtual bool mightVanish() const = 0;
};

/** \brief Wrapper of a Dune entity of type DuneEntity and codimension codim

 \note The codimension must be given explicitly (even though it could be
 derived from the traits of DuneEntity) because this class needs to be
 specialised for entities of codimension 0.
 */
template<int codim, typename DuneEntity>
class ConcreteEntity: public Entity<codim> {
	dune_static_assert((int)DuneEntity::codimension == (int)Entity<codim>::codimension,
			"ConcreteEntity: codimension mismatch");

protected:
	const DuneEntity* m_dune_entity;
	/** \internal Entity geometry. Updated on demand (on calling
	 * geometry()), hence declared as mutable. */
	mutable ConcreteGeometry<typename DuneEntity::Geometry> m_geometry;

	template<typename > friend class ConcreteEntityPointer;
	template<typename > friend class ConcreteRangeEntityIterator;
	template<typename, int> friend class ConcreteSubentityIterator;

	void setDuneEntity(const DuneEntity* dune_entity) {
		m_dune_entity = dune_entity;
	}

public:
	/** Default constructor */
	ConcreteEntity() :
			m_dune_entity(0) {
	}

	/** Constructor from a pointer to DuneEntity */
	ConcreteEntity(const DuneEntity* dune_entity) :
			m_dune_entity(dune_entity) {
	}

	const DuneEntity& duneEntity() const {
		return *m_dune_entity;
	}

	virtual int level() const {
		return m_dune_entity->level();
	}

	virtual const Geometry& geometry() const {
		m_geometry.setDuneGeometry(&m_dune_entity->geometry());
		return m_geometry;
	}

	virtual Dune::GeometryType type() const {
		return m_dune_entity->type();
	}
};

/** \brief Wrapper of a Dune entity of type DuneEntity and codimension 0
 */

template<typename DuneEntity>
class ConcreteEntity<0, DuneEntity> : public Entity<0> {
	dune_static_assert((int)DuneEntity::codimension == (int)codimension,
			"ConcreteEntity: codimension mismatch");

protected:
	const DuneEntity* m_dune_entity;
	/** \internal Entity geometry. Updated on demand (on calling
	 * geometry()), hence declared as mutable. */
	mutable ConcreteGeometry<typename DuneEntity::Geometry> m_geometry;

	template<typename > friend class ConcreteEntityPointer;
	template<typename > friend class ConcreteRangeEntityIterator;
	template<typename, int> friend class ConcreteSubentityIterator;

	void setDuneEntity(const DuneEntity* dune_entity) {
		m_dune_entity = dune_entity;
	}

public:
	/** Default constructor */
	ConcreteEntity() :
			m_dune_entity(0) {
	}

	/** Constructor from a pointer to DuneEntity */
	explicit ConcreteEntity(const DuneEntity* dune_entity) :
			m_dune_entity(dune_entity) {
	}

	/** Read-only access to the underlying DuneEntity object */
	const DuneEntity& duneEntity() const {
		return *m_dune_entity;
	}

	virtual int level() const {
		return m_dune_entity->level();
	}

	virtual const Geometry& geometry() const {
		m_geometry.setDuneGeometry(&m_dune_entity->geometry());
		return m_geometry;
	}

	virtual Dune::GeometryType type() const {
		return m_dune_entity->type();
	}

	virtual int edgeCount() const {
		return m_dune_entity->template count<1>();
	}

	virtual int vertexCount() const {
		return m_dune_entity->template count<2>();
	}

	virtual EntityIterator<1>* subEdgeIterator() const;
	virtual EntityIterator<2>* subVertexIterator() const;

// TODO: write properly these definitions
//  virtual SurfaceGridIntersectionIterator ileafbegin() const
//  {
//    typedef typename GridType::LeafIntersectionIterator DuneIterator;
//    typedef ConcreteSurfaceGridIntersectionIterator<DuneIterator> ConcreteIterator;
//    DuneIterator dit = m_dune_entity->ileafbegin();
//    return SurfaceGridIntersectionIterator(new ConcreteIterator(dit));
//  }

//  virtual SurfaceGridIntersectionIterator ileafend() const
//  {
//    typedef typename GridType::LeafIntersectionIterator DuneIterator;
//    typedef ConcreteSurfaceGridIntersectionIterator<DuneIterator> ConcreteIterator;
//    DuneIterator dit = m_dune_entity->ileafend();
//    return SurfaceGridIntersectionIterator(new ConcreteIterator(dit));
//  }

//  virtual SurfaceGridIntersectionIterator ilevelbegin() const
//  {
//    typedef typename GridType::LevelIntersectionIterator DuneIterator;
//    typedef ConcreteSurfaceGridIntersectionIterator<DuneIterator> ConcreteIterator;
//    DuneIterator dit = m_dune_entity->ilevelbegin();
//    return SurfaceGridIntersectionIterator(new ConcreteIterator(dit));
//  }

//  virtual SurfaceGridIntersectionIterator ilevelend() const
//  {
//    typedef typename GridType::LevelIntersectionIterator DuneIterator;
//    typedef ConcreteSurfaceGridIntersectionIterator<DuneIterator> ConcreteIterator;
//    DuneIterator dit = m_dune_entity->ilevelend();
//    return SurfaceGridIntersectionIterator(new ConcreteIterator(dit));
//  }

	virtual EntityPointer<0>* father() const;

	virtual bool hasFather() const {
		return m_dune_entity->hasFather();
	}

	virtual bool isLeaf() const {
		return m_dune_entity->isLeaf();
	}

	virtual bool isRegular() const {
		return m_dune_entity->isRegular();
	}

// To be implemented.
//  virtual const SurfaceGridLocalGeometry& geometryInFather() const
//  { return SurfaceGridLocalGeometry(m_dune_entity->geometryInFather()); }

	virtual EntityIterator<0>* sonIterator(int maxlevel) const;

	virtual bool isNew() const {
		return m_dune_entity->isNew();
	}

	virtual bool mightVanish() const {
		return m_dune_entity->mightVanish();
	}
};

} // namespace Bempp

#endif
