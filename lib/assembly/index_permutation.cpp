#include "index_permutation.hpp"

#include "../common/boost_make_shared_fwd.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_LocalMap.h>
#include <Epetra_SerialComm.h>

namespace Bempp {

shared_ptr<const Epetra_CrsMatrix> IndexPermutation::permutationMatrix() const {
  const int size = m_permutedIndices.size();
  Epetra_SerialComm comm; // To be replaced once we begin to use MPI
  Epetra_LocalMap rowMap(size, 0 /* index_base */, comm);
  Epetra_LocalMap colMap(size, 0 /* index_base */, comm);
  shared_ptr<Epetra_CrsMatrix> result = boost::make_shared<Epetra_CrsMatrix>(
      Copy, rowMap, colMap, 1 /* one entry per row */);

  const double ONE = 1.;
  for (int i = 0; i < size; ++i)
    result->InsertGlobalValues(m_permutedIndices[i], 1 /* one entry */, &ONE,
                               &i);
  result->FillComplete();
  return result;
}

} // namespace Bempp
