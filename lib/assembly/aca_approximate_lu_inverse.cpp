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

#include "bempp/common/config_ahmed.hpp"

#ifdef WITH_AHMED
#include "aca_approximate_lu_inverse.hpp"

#include "ahmed_aux.hpp"
#include "discrete_aca_boundary_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/serial_blas_region.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_DetachedSpmdVectorView.hpp>
#include <Thyra_DefaultSpmdVectorSpace_decl.hpp>
#endif

#include <tbb/tick_count.h>

namespace Bempp {

template <typename ValueType>
AcaApproximateLuInverse<ValueType>::AcaApproximateLuInverse(
    const DiscreteAcaBoundaryOperator<ValueType> &fwdOp, MagnitudeType delta,
    VerbosityLevel::Level verbosityLevel)
    :
// All range-domain swaps intended!
#ifdef WITH_TRILINOS
      m_domainSpace(fwdOp.m_rangeSpace),
      m_rangeSpace(fwdOp.m_domainSpace),
#else
      m_rowCount(fwdOp.columnCount()),
      m_columnCount(fwdOp.rowCount()),
#endif
      m_blockCluster(0), m_blocksL(0), m_blocksU(0),
      m_domainPermutation(fwdOp.m_rangePermutation),
      m_rangePermutation(fwdOp.m_domainPermutation) {
  const bool verbosityAtLeastDefault =
      (verbosityLevel >= VerbosityLevel::DEFAULT);
  if (verbosityAtLeastDefault)
    std::cout << "Starting H-LU decomposition..." << std::endl;
  tbb::tick_count start = tbb::tick_count::now();
  const blcluster *fwdBlockCluster = fwdOp.m_blockCluster.get();
  bool result = genLUprecond(const_cast<blcluster *>(fwdBlockCluster),
                             fwdOp.m_blocks.get(), delta, fwdOp.m_maximumRank,
                             m_blockCluster, m_blocksL, m_blocksU, true);
  tbb::tick_count end = tbb::tick_count::now();
  if (!result)
    throw std::runtime_error(
        "AcaApproximateLuInverse::AcaApproximateLuInverse(): "
        "Approximate LU factorisation failed");

  if (verbosityAtLeastDefault) {
    std::cout << "H-LU decomposition took " << (end - start).seconds() << " s"
              << std::endl;
    size_t origMemory =
        sizeof(ValueType) * m_domainSpace->dim() * m_rangeSpace->dim();
    size_t ahmedMemory = sizeH(m_blockCluster, m_blocksL, 'L') +
                         sizeH(m_blockCluster, m_blocksU, 'U');
    int maximumRankL = Hmax_rank(m_blockCluster, m_blocksL, 'L');
    int maximumRankU = Hmax_rank(m_blockCluster, m_blocksU, 'U');
    std::cout << "\nNeeded storage: " << ahmedMemory / 1024. / 1024. << " MB.\n"
              << "Without approximation: " << origMemory / 1024. / 1024.
              << " MB.\n"
              << "Compressed to " << (100. * ahmedMemory) / origMemory << "%.\n"
              << "Maximum rank: " << std::max(maximumRankL, maximumRankU)
              << ".\n" << std::endl;
  }
}

template <>
AcaApproximateLuInverse<std::complex<float>>::AcaApproximateLuInverse(
    const DiscreteAcaBoundaryOperator<std::complex<float>> &fwdOp,
    MagnitudeType delta, VerbosityLevel::Level verbosityLevel)
    :
// All range-domain swaps intended!
#ifdef WITH_TRILINOS
      m_domainSpace(fwdOp.m_rangeSpace),
      m_rangeSpace(fwdOp.m_domainSpace),
#else
      m_rowCount(fwdOp.columnCount()),
      m_columnCount(fwdOp.rowCount()),
#endif
      m_blockCluster(0), m_blocksL(0), m_blocksU(0),
      m_domainPermutation(fwdOp.m_rangePermutation),
      m_rangePermutation(fwdOp.m_domainPermutation) {
  // Ahmed doesn't define the genLUprecond() variant with
  // the second parameter of type mblock<scomp>**
  throw std::runtime_error(
      "AcaApproximateLuInverse::AcaApproximateLuInverse(): "
      "due to a deficiency in Ahmed approximate LU factorisation "
      "of single-precision complex H matrices is not supported");
}

template <typename ValueType>
AcaApproximateLuInverse<ValueType>::~AcaApproximateLuInverse() {
  if (m_blockCluster) {
    freembls(m_blockCluster, m_blocksL);
    freembls(m_blockCluster, m_blocksU);
    delete m_blockCluster;
  }
}

template <typename ValueType>
unsigned int AcaApproximateLuInverse<ValueType>::rowCount() const {
#ifdef WITH_TRILINOS
  return m_rangeSpace->dim();
#else
  return m_rowCount;
#endif
}

template <typename ValueType>
unsigned int AcaApproximateLuInverse<ValueType>::columnCount() const {
#ifdef WITH_TRILINOS
  return m_domainSpace->dim();
#else
  return m_columnCount;
#endif
}

template <typename ValueType>
void AcaApproximateLuInverse<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, arma::Mat<ValueType> &block) const {
  throw std::runtime_error("AcaApproximateLuInverse::addBlock(): "
                           "not implemented");
}

#ifdef WITH_TRILINOS

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
AcaApproximateLuInverse<ValueType>::domain() const {
  return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
AcaApproximateLuInverse<ValueType>::range() const {
  return m_rangeSpace;
}

template <typename ValueType>
bool AcaApproximateLuInverse<ValueType>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  // TODO: implement remaining variants (transpose & conjugate transpose)
  return (M_trans == Thyra::NOTRANS);
}
#endif // WITH_TRILINOS

template <typename ValueType>
void AcaApproximateLuInverse<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const arma::Col<ValueType> &x_in,
    arma::Col<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {
  if (trans != NO_TRANSPOSE)
    throw std::runtime_error(
        "AcaApproximateLuInverse::applyBuiltInImpl(): "
        "transposition modes other than NO_TRANSPOSE are not supported");
  if (columnCount() != x_in.rows() && rowCount() != y_inout.rows())
    throw std::invalid_argument("AcaApproximateLuInverse::applyBuiltInImpl(): "
                                "incorrect vector length");

  if (beta == static_cast<ValueType>(0.))
    y_inout.fill(static_cast<ValueType>(0.));
  else
    y_inout *= beta;

  // will act both as a permuted argument and permuted result
  arma::Col<ValueType> permuted;
  m_domainPermutation.permuteVector(x_in, permuted);

  HLU_solve(m_blockCluster, m_blocksL, m_blocksU, ahmedCast(permuted.memptr()));

  arma::Col<ValueType> operatorActionResult;
  m_rangePermutation.unpermuteVector(permuted, operatorActionResult);
  y_inout += alpha * operatorActionResult;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(AcaApproximateLuInverse);

} // namespace Bempp

#endif // WITH_AHMED
