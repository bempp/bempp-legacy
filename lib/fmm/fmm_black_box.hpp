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

#ifndef bempp_black_box_hpp
#define bempp_black_box_hpp

#include <vector>

#include "../common/armadillo_fwd.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"


namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename KernelType> class CollectionOfKernels;
template <typename KernelFunctor> class DefaultCollectionOfKernels;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename ValueType> class FmmTransform;
/** \endcond */

template <typename KernelType, typename ValueType>
class FmmBlackBox : public FmmTransform<ValueType>
{
public:
	typedef typename FmmTransform<ValueType>::CoordinateType CoordinateType;
	typedef Fiber::CollectionOfKernels<KernelType> CollectionOfKernels;

	template <typename KernelFunctor>
	FmmBlackBox(const KernelFunctor& kernelFunctor, unsigned int n);

	// multipole to multipole (M2M) translation matrix
	virtual arma::Mat<ValueType> M2M(
		const arma::Col<CoordinateType>& childPosition, 
		const arma::Col<CoordinateType>& parentPosition) const;

	// multipole to local (M2L) translation matrix
	virtual arma::Mat<ValueType> M2L(
		const arma::Col<CoordinateType>& sourceCentre, 
		const arma::Col<CoordinateType>& fieldCentre,
		const arma::Col<CoordinateType>& boxSize) const;

	// local to local (L2L) translation matrix
	virtual arma::Mat<ValueType> L2L(
		const arma::Col<CoordinateType>& parentPosition, 
		const arma::Col<CoordinateType>& childPosition) const;

	virtual void generateGaussPoints();

	virtual void getKernelWeight(
		arma::Mat<ValueType>& kernelWeightMat,
		arma::Col<ValueType>& kernelWeightVec) const;

	unsigned int getN() const {return m_n;}

protected:

	virtual void evaluateAtGaussPointS(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const;

	virtual void evaluateAtGaussPointDiffS(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const;

private:
	unsigned int m_n;
	arma::Mat<CoordinateType> m_Tk;
	shared_ptr<CollectionOfKernels> m_kernels;
};

// Constructor cannot be placed in fmm_black_box.hpp, otherwise a linker error 
// results based on general_elementary_singular_integral_operator_imp.hpp. Note 
// that the very reason for the presence of the imp files is to avoid linker 
// errors. They are (and must be) directly included, yet they contain code.
template <typename KernelType, typename ValueType>
template <typename KernelFunctor>
FmmBlackBox<KernelType, ValueType>::FmmBlackBox(
	const KernelFunctor& kernelFunctor, unsigned int n)
 :	m_kernels(
        new Fiber::DefaultCollectionOfKernels<KernelFunctor>(kernelFunctor)),
	m_n(n), m_Tk(n, n), FmmTransform<ValueType>(n*n*n, true)
{
	generateGaussPoints();
}

template <typename KernelType, typename ValueType>
class FmmSingleLayerBlackBox : public FmmBlackBox<KernelType, ValueType>
{
public:
	typedef typename FmmBlackBox<KernelType, ValueType>::CoordinateType CoordinateType;

	template <typename KernelFunctor>
	FmmSingleLayerBlackBox(const KernelFunctor& kernelFunctor, unsigned int n)
		: FmmBlackBox<KernelType, ValueType>(kernelFunctor, n) {}

	virtual void evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const;

	virtual void evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const;
};

template <typename KernelType, typename ValueType>
class FmmDoubleLayerBlackBox : public FmmBlackBox<KernelType, ValueType>
{
public:
	typedef typename FmmBlackBox<KernelType, ValueType>::CoordinateType CoordinateType;

	template <typename KernelFunctor>
	FmmDoubleLayerBlackBox(const KernelFunctor& kernelFunctor, unsigned int n)
		: FmmBlackBox<KernelType, ValueType>(kernelFunctor, n) {}

	virtual void evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const;

	virtual void evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const;
};

template <typename KernelType, typename ValueType>
class FmmAdjointDoubleLayerBlackBox : public FmmBlackBox<KernelType, ValueType>
{
public:
	typedef typename FmmBlackBox<KernelType, ValueType>::CoordinateType CoordinateType;

	template <typename KernelFunctor>
	FmmAdjointDoubleLayerBlackBox(const KernelFunctor& kernelFunctor, unsigned int n)
		: FmmBlackBox<KernelType, ValueType>(kernelFunctor, n) {}

	virtual void evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const;

	virtual void evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const;
};

template <typename KernelType, typename ValueType>
class FmmHypersingularBlackBox : public FmmBlackBox<KernelType, ValueType>
{
public:
	typedef typename FmmBlackBox<KernelType, ValueType>::CoordinateType CoordinateType;

	template <typename KernelFunctor>
	FmmHypersingularBlackBox(const KernelFunctor& kernelFunctor, unsigned int n)
		: FmmBlackBox<KernelType, ValueType>(kernelFunctor, n) {}

	virtual void evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const;

	virtual void evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const;
};


} // namespace Bempp

#endif
