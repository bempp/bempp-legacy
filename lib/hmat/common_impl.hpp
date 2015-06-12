// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_COMMON_IMPL_HPP
#define HMAT_COMMON_IMPL_HPP

#include "common.hpp"

namespace hmat {

inline IndexSetType fillIndexRange(std::size_t start, std::size_t stop) {
  IndexSetType result;
  result.reserve(stop - start);
  for (int i = start; i < stop; ++i)
    result.push_back(i);
  return result;
}

}

#endif
