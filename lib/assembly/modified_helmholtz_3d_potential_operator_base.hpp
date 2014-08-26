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

#ifndef bempp_modified_helmholtz_3d_potential_operator_base_hpp
#define bempp_modified_helmholtz_3d_potential_operator_base_hpp

#include "elementary_potential_operator.hpp"

#include <boost/scoped_ptr.hpp>

namespace Bempp {

/** \ingroup helmholtz_3d
 *  \brief Base class for potential operators for the modified Helmholtz
 *equation in 3D.
 *
 *  \tparam Impl
 *    Type of the internal implementation object.
 *  \tparam BasisFunctionType_
 *    Type of the values of the basis functions into which functions acted upon
 *    by the operator are expanded. It can take the following values: \c float,
 *    \c double, <tt>std::complex<float></tt> and
 *    <tt>std::complex<double></tt>.
 *
 *  \see helmholtz_3d */
template <typename Impl, typename BasisFunctionType_>
class ModifiedHelmholtz3dPotentialOperatorBase
    : public ElementaryPotentialOperator<
          BasisFunctionType_,
          typename ScalarTraits<BasisFunctionType_>::ComplexType,
          typename ScalarTraits<BasisFunctionType_>::ComplexType> {
  typedef ElementaryPotentialOperator<
      BasisFunctionType_,
      typename ScalarTraits<BasisFunctionType_>::ComplexType,
      typename ScalarTraits<BasisFunctionType_>::ComplexType> Base;

public:
  /** \brief Type of the values of the basis functions into which functions
   *  acted upon by the operator are expanded. */
  typedef typename Base::BasisFunctionType BasisFunctionType;
  /** \brief Type of the values of kernel functions. */
  typedef typename Base::KernelType KernelType;
  /** \brief Type of the values of the potential. */
  typedef typename Base::ResultType ResultType;
  /** \copydoc ElementaryPotentialOperator::CoordinateType */
  typedef typename Base::CoordinateType CoordinateType;
  /** \copydoc ElementaryPotentialOperator::CollectionOfShapesetTransformations
   */
  typedef typename Base::CollectionOfShapesetTransformations
  CollectionOfShapesetTransformations;
  /** \copydoc ElementaryPotentialOperator::CollectionOfBasisTransformations */
  typedef typename Base::CollectionOfBasisTransformations
  CollectionOfBasisTransformations;
  /** \copydoc ElementaryPotentialOperator::CollectionOfKernels */
  typedef typename Base::CollectionOfKernels CollectionOfKernels;
  /** \copydoc ElementaryPotentialOperator::KernelTrialIntegral */
  typedef typename Base::KernelTrialIntegral KernelTrialIntegral;

  /** \brief Constructor.
   *
   *  \param[in] waveNumber
   *    Wave number. See \ref modified_helmholtz_3d for its definition. */
  ModifiedHelmholtz3dPotentialOperatorBase(KernelType waveNumber);
  /** \brief Copy constructor. */
  ModifiedHelmholtz3dPotentialOperatorBase(
      const ModifiedHelmholtz3dPotentialOperatorBase &other);
  /** \brief Destructor. */
  virtual ~ModifiedHelmholtz3dPotentialOperatorBase();
  /** \brief Assignment operator. */
  ModifiedHelmholtz3dPotentialOperatorBase &
  operator=(const ModifiedHelmholtz3dPotentialOperatorBase &rhs);

  /** \brief Return the wave number set previously in the constructor. */
  KernelType waveNumber() const;

private:
  virtual const CollectionOfKernels &kernels() const;
  virtual const CollectionOfShapesetTransformations &
  trialTransformations() const;
  virtual const KernelTrialIntegral &integral() const;

private:
  /** \cond PRIVATE */
  boost::scoped_ptr<Impl> m_impl;
  /** \endcond */
};

} // namespace Bempp

#endif
