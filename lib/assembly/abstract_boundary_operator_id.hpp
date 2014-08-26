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

#ifndef bempp_abstract_boundary_operator_id_hpp
#define bempp_abstract_boundary_operator_id_hpp

#include "../common/common.hpp"
#include "../common/deprecated.hpp"
#include "../common/shared_ptr.hpp"

#include <iostream>
#include <boost/functional/hash.hpp>

namespace Bempp {

/** \ingroup abstract_boundary_operators
 *  \brief Base class for identifiers of an abstract boundary operator.
 *
 *  \deprecated This class is deprecated and will be removed in a
 *  future version of BEM++. Boundary operator identifiers are no longer used
 *  by the discrete weak-form caching mechanism.
 */
class AbstractBoundaryOperatorId {
public:
  virtual ~AbstractBoundaryOperatorId() {}
  virtual size_t hash() const = 0;
  virtual void dump() const {}
  virtual bool isEqual(const AbstractBoundaryOperatorId &other) const = 0;

  bool operator==(const AbstractBoundaryOperatorId &other) const {
    return isEqual(other);
  }

  bool operator!=(const AbstractBoundaryOperatorId &other) const {
    return !operator==(other);
  }
};

} // namespace Bempp

namespace tbb {

inline size_t tbb_hasher(
    const boost::shared_ptr<const Bempp::AbstractBoundaryOperatorId> &id) {
  return id->hash();
}

inline size_t tbb_hasher(float id) { return boost::hash<float>()(id); }

inline size_t tbb_hasher(double id) { return boost::hash<double>()(id); }

} // namespace tbb

#include <tbb/concurrent_unordered_map.h>

namespace Bempp {

// copied from Boost
template <typename T> inline void tbb_hash_combine(size_t &seed, const T &v) {
  seed ^= tbb::tbb_hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
}

#endif // bempp_abstract_boundary_operator_id_hpp
