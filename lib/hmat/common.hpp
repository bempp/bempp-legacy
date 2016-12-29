// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_COMMON_HPP
#define HMAT_COMMON_HPP

#include "shared_ptr.hpp"

#include "point.hpp"
#include <array>
#include <vector>

#include <mpi.h>

namespace hmat {

typedef std::array<std::size_t, 4> BlockIndexRangeType;
typedef std::array<std::size_t, 2> IndexRangeType;
typedef std::vector<std::size_t> IndexSetType;

enum RowColSelector { ROW, COL };

enum TransposeMode { NOTRANS, TRANS, CONJ, CONJTRANS };

enum DataBlockType {

  DENSE,
  LOW_RANK_AB

};

IndexSetType fillIndexRange(std::size_t start, std::size_t stop);
}

// MPI Types

template <typename T> struct MpiTrait {};

template <> struct MpiTrait<int> {

  inline MpiTrait() : type(MPI_INT) {}

  MPI_Datatype type;
};

template <> struct MpiTrait<float> {

  inline MpiTrait() : type(MPI_FLOAT) {}

  MPI_Datatype type;
};

template <> struct MpiTrait<double> {

  inline MpiTrait() : type(MPI_DOUBLE) {}

  MPI_Datatype type;
};

template <> struct MpiTrait<std::complex<float>> {

  inline MpiTrait() : type(MPI_C_FLOAT_COMPLEX) {}

  MPI_Datatype type;
};

template <> struct MpiTrait<std::complex<double>> {

  inline MpiTrait() : type(MPI_C_DOUBLE_COMPLEX) {}

  MPI_Datatype type;
};

#include "common_impl.hpp"

#endif
