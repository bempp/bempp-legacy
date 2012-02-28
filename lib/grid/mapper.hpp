#ifndef bempp_mapper_hpp
#define bempp_mapper_hpp

namespace Bempp
{

template <int codim> class Entity;

/** \brief Abstract mapper class.

  A mapper provides a mapping from a set of entities to a range of integers from
  0 to (size() - 1), where size() is the number of entities contained in the set.
  */
class Mapper
{
public:    
    /** \brief Total number of entities in the entity set managed by the mapper. */
    virtual int size() const = 0;

    /** \brief Index of the entity \e of codimension 0.

     The result of calling this method with an entity that does not belong to
     the mapped set is undefined.

     \return An index in the range 0 ... (size() - 1). */
    virtual int entityIndex(const Entity<0>& e) const = 0;
    /** \brief Index of the entity \e of codimension 1.

     \overload
     */
    virtual int entityIndex(const Entity<1>& e) const = 0;
    /** \brief Index of the entity \e of codimension 2.

     \overload
     */
    virtual int entityIndex(const Entity<2>& e) const = 0;
    /** \brief Index of the entity \e of codimension 3.

     \overload
     */
    virtual int entityIndex(const Entity<3>& e) const = 0;
    /** \brief Index of \p i'th subentity of codimension \p codimSub of
     entity \p e of codimension 0. */

    virtual int subEntityIndex(const Entity<0>& e, int i,
                               unsigned int codimSub) const = 0;
};

} // namespace Bempp

#endif
