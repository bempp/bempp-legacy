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

#ifndef bempp_concrete_grid_view_hpp
#define bempp_concrete_grid_view_hpp

#include "grid_view.hpp"
#include "concrete_entity.hpp"
#include "concrete_index_set.hpp"
#include "concrete_range_entity_iterator.hpp"
#include "concrete_vtk_writer.hpp"
#include "geometry_type.hpp"

namespace Bempp
{

/** \brief Wrapper of a Dune grid view of type \p DuneGridView. */
template<typename DuneGridView>
class ConcreteGridView: public GridView
{
private:
    DuneGridView m_dune_gv;
    ConcreteIndexSet<DuneGridView> m_index_set;

public:
    /** \brief Constructor */
    explicit ConcreteGridView(const DuneGridView& dune_gv) :
        m_dune_gv(dune_gv), m_index_set(&dune_gv.indexSet()) {
    }

    /** \brief Read-only access to the underlying Dune grid view object. */
    const DuneGridView& duneGridView() const {
        return m_dune_gv;
    }

    /** \brief Access to the underlying Dune grid view object. Use at your own risk! */
    DuneGridView& duneGridView() {
        return m_dune_gv;
    }

    virtual const IndexSet& indexSet() const {
        return m_index_set;
    }

    virtual int entityCount(int codim) const {
        return m_dune_gv.size(codim);
    }

    virtual int entityCount(const GeometryType &type) const {
        return m_dune_gv.size(type);
    }

    virtual bool containsEntity(const Entity<0>& e) const {
        return containsEntityCodimN(e);
    }
    virtual bool containsEntity(const Entity<1>& e) const {
        return containsEntityCodimN(e);
    }
    virtual bool containsEntity(const Entity<2>& e) const {
        return containsEntityCodimN(e);
    }
    virtual bool containsEntity(const Entity<3>& e) const {
        return containsEntityCodimN(e);
    }

    virtual std::auto_ptr<VtkWriter> vtkWriter(Dune::VTK::DataMode dm=Dune::VTK::conforming) const {
        return std::auto_ptr<VtkWriter>(new ConcreteVtkWriter<DuneGridView>(m_dune_gv, dm));
    }

private:
    virtual std::auto_ptr<EntityIterator<0> > entityCodim0Iterator() const {
        return entityCodimNIterator<0>();
    }
    virtual std::auto_ptr<EntityIterator<1> > entityCodim1Iterator() const {
        return entityCodimNIterator<1>();
    }
    virtual std::auto_ptr<EntityIterator<2> > entityCodim2Iterator() const {
        return entityCodimNIterator<2>();
    }
    virtual std::auto_ptr<EntityIterator<3> > entityCodim3Iterator() const {
        return entityCodimNIterator<3>();
    }

    template <int codim>
    typename boost::disable_if_c<codim <= DuneGridView::dimension, bool>::type
    containsEntityCodimN(const Entity<codim>& e) const {
        throw std::logic_error("GridView::containsEntity(): invalid entity codimension");
    }

    template <int codim>
    typename boost::enable_if_c<codim <= DuneGridView::dimension, bool>::type
    containsEntityCodimN(const Entity<codim>& e) const {
        typedef typename DuneGridView::template Codim<codim>::Entity DuneEntity;
        typedef ConcreteEntity<codim, DuneEntity> ConcEntity;
        const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
        return m_dune_gv.contains(ce.duneEntity());
    }

    template <int codim>
    typename boost::disable_if_c<codim <= DuneGridView::dimension, std::auto_ptr<EntityIterator<codim> > >::type
    entityCodimNIterator() const {
        throw std::logic_error("GridView::entityIterator(): invalid entity codimension");
    }

    template <int codim>
    typename boost::enable_if_c<codim <= DuneGridView::dimension, std::auto_ptr<EntityIterator<codim> > >::type
    entityCodimNIterator() const {
        typedef typename DuneGridView::template Codim<codim>::Iterator DuneIterator;
        typedef typename DuneGridView::template Codim<codim>::EntityPointer DuneEntityPointer;
        typedef ConcreteRangeEntityIterator<DuneIterator, DuneEntityPointer> ConcIterator;
        return std::auto_ptr<EntityIterator<codim> >(
                   new ConcIterator(m_dune_gv.template begin<codim>(),
                                    m_dune_gv.template end<codim>()));
    }
};

} // namespace Bempp

#endif
