// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_SHARED_PTR_HPP
#define HMAT_SHARED_PTR_HPP

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <iostream>

#include <tbb/concurrent_unordered_map.h>

namespace hmat {

using boost::shared_ptr;
using boost::make_shared;
using boost::enable_shared_from_this;
using boost::weak_ptr;
using boost::const_pointer_cast;
using boost::static_pointer_cast;

template <typename T> struct shared_ptr_hash {

  std::size_t operator()(const shared_ptr<T> &key) const {
    return tbb::tbb_hasher(key.get());
  }
};
}

#endif
