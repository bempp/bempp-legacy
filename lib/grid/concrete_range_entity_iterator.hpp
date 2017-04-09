// Copyright (C) 2011-2012 by the BEM++ Authors
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

#ifndef bempp_concrete_range_entity_iterator_hpp
#define bempp_concrete_range_entity_iterator_hpp

#include "../common/common.hpp"

#include "entity_iterator.hpp"
#include "concrete_entity_decl.hpp"

namespace Bempp {

class DomainIndex;

/** \ingroup grid_internal
\brief Iterator over entities referenced by a range of Dune iterators of
type \p DuneEntityIt.

\p DuneEntityPointer is the base class of \p DuneEntityIt (it would be easily
read from the DuneEntityIt::Base typedef, but, alas, this typedef is
private...).
*/
template <typename DuneEntityIt>
class ConcreteRangeEntityIterator
    : public EntityIterator<DuneEntityIt::codimension> {
private:
  DuneEntityIt m_begin, m_end, m_cur;
  std::unique_ptr<ConcreteEntity<ConcreteRangeEntityIterator::codimension,
                 typename DuneEntityIt::Entity>> m_entity;
  const DomainIndex &m_domain_index;

  void updateEntity() {
    if (!this->finished())
      m_entity.reset(new ConcreteEntity<ConcreteRangeEntityIterator::codimension,
                                 typename DuneEntityIt::Entity>(*m_cur, m_domain_index));
  }

  void updateFinished() { this->m_finished = (m_cur == m_end); }

public:
  /** \brief Constructor. The iterator will go over the range [\p begin, \p
   * end). */
  ConcreteRangeEntityIterator(const DuneEntityIt &begin,
                              const DuneEntityIt &end,
                              const DomainIndex &domain_index)
      : m_begin(begin), m_end(end), m_cur(begin),
        m_entity(new ConcreteEntity<ConcreteRangeEntityIterator::codimension,
                 typename DuneEntityIt::Entity>(*begin, domain_index)),
        m_domain_index(domain_index) {
    updateFinished();
  }

  virtual void next() {
    ++m_cur;
    updateFinished();
    updateEntity();
  }

  virtual const Entity<ConcreteRangeEntityIterator::codimension> &
  entity() const {
    return *m_entity;
  }

  virtual std::unique_ptr<
      EntityPointer<ConcreteRangeEntityIterator::codimension>>
  frozen() const {
    const int codim = ConcreteRangeEntityIterator::codimension;
    return std::unique_ptr<EntityPointer<codim>>(
        new ConcreteEntityPointer<typename DuneEntityIt::Entity>(*m_cur, m_domain_index));
  }
};

} // namespace Bempp

#endif
