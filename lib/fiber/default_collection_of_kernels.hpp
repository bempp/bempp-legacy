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

#ifndef fiber_default_collection_of_kernels_hpp
#define fiber_default_collection_of_kernels_hpp

#include "collection_of_kernels.hpp"

namespace Fiber {

/** \ingroup weak_form_elements
 *  \brief Default implementation of a collection of kernels.

    This class implements the interface defined by CollectionOfKernels using
    a functor object to evaluate the kernels at specific point pairs.

    \tparam Functor
      Type of the functor that will be passed to the constructor and used to
      evaluate a number of kernels at individual point pairs.

    The Functor class should provide the following interface:

    \code
    class KernelCollectionFunctor
    {
    public:
        typedef ... ValueType; // can be float, double, std::complex<float>
                               // or std::complex<double>
        typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

        // Return the number of kernels in the collection
        int kernelCount() const;
        // Return the number of rows in the tensor being the value of i'th
 kernel
        int kernelRowCount(int i) const;
        // Return the number of columns in the tensor being the value of i'th
 kernel
        int kernelColCount(int i) const;

        // Specify types of geometrical data required by the kernels (see
        // documentation of CollectionOfKernels::addGeometricalDependencies
        // for details)
        void addGeometricalDependencies(size_t& testGeomDeps,
                                        size_t& trialGeomDeps) const;

        // Evaluate the kernels at the pair of points whose geometrical data
        // are provided in the testGeomData and trialGeomData arguments. The
        // (j, k)th element of the tensor being the value of i'th kernel should
        // be written to result[i](j, k).
        template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
        void evaluate(
                const ConstGeometricalDataSlice<CoordinateType>& testGeomData,
                const ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
                CollectionOf2dSlicesOfNdArrays<ValueType>& result) const;

        // (Optional)
        // Return an estimate of the magnitude of the kernel at test and trial
        // points lying in a given distance from each other. This estimate does
        // not need to be accurate; in practice, it is only useful to give it at
        // all for exponentially decaying kernels.  If this function is not
        // defined, the kernel behaves as if its estimated magnitude was 1
        // everywhere.
        CoordinateType estimateRelativeScale(CoordinateType distance) const;
    };
    \endcode

    See the Laplace3dSingleLayerPotentialKernelFunctor class for an example
    implementation of a (simple) kernel collection functor.
 */
template <typename Functor>
class DefaultCollectionOfKernels
    : public CollectionOfKernels<typename Functor::ValueType> {
  typedef CollectionOfKernels<typename Functor::ValueType> Base;

public:
  typedef typename Base::ValueType ValueType;
  typedef typename Base::CoordinateType CoordinateType;

  explicit DefaultCollectionOfKernels(const Functor &functor)
      : m_functor(functor) {}

  const Functor &functor() const { return m_functor; }

  Functor &functor() { return m_functor; }

  virtual void addGeometricalDependencies(size_t &testGeomDeps,
                                          size_t &trialGeomDeps) const;

  virtual void
  evaluateAtPointPairs(const GeometricalData<CoordinateType> &testGeomData,
                       const GeometricalData<CoordinateType> &trialGeomData,
                       CollectionOf3dArrays<ValueType> &result) const;

  virtual void
  evaluateOnGrid(const GeometricalData<CoordinateType> &testGeomData,
                 const GeometricalData<CoordinateType> &trialGeomData,
                 CollectionOf4dArrays<ValueType> &result) const;

  virtual std::pair<const char *, int> evaluateClCode() const;

  virtual CoordinateType estimateRelativeScale(CoordinateType distance) const;

private:
  Functor m_functor;
};

} // namespace Fiber

#include "default_collection_of_kernels_imp.hpp"

#endif
