// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_COMMON_HPP
#define HMAT_COMMON_HPP

#include "shared_ptr.hpp"

#include <array>
#include <vector>

namespace hmat {

typedef std::array<std::size_t, 4> BlockIndexRangeType;
typedef std::array<std::size_t, 2> IndexRangeType;
typedef std::vector<std::size_t> IndexSetType;

enum class TransposeMode {
  NOTRANS,
  TRANS,
  CONJTRANS
};

IndexSetType fillIndexRange(std::size_t start, std::size_t stop);
}

#include "common_impl.hpp"

#endif
