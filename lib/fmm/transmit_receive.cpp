
#include "transmit_receive.hpp"
#include "fmm_transform.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../assembly/surface_normal_independent_function.hpp"
#include "../assembly/surface_normal_dependent_function.hpp"

namespace Bempp
{

// see examples/tutorial_dirichlet.cpp
// and fiber/surface_normal_independent_function.hpp
// returns the TransmitFunction (Antenna) function for a given Gauss point khat

// trial function

template <typename ResultType>
FmmFunctionMultiplyingTrial<ResultType>::FmmFunctionMultiplyingTrial(
		const arma::Col<CoordinateType>& khat, 
		const arma::Col<CoordinateType>& nodeCentre,
		const arma::Col<CoordinateType>& nodeSize,
		const FmmTransform<ResultType>& fmmTransform)
	 :	m_khat(khat), m_nodeCentre(nodeCentre), m_nodeSize(nodeSize), 
		m_fmmTransform(fmmTransform)
{
}

template <typename ResultType>
int FmmFunctionMultiplyingTrial<ResultType>::argumentDimension() const
{
	return 3;
}

template <typename ResultType>
int FmmFunctionMultiplyingTrial<ResultType>::resultDimension() const 
{
	return 1;
}

template <typename ResultType>
inline void FmmFunctionMultiplyingTrial<ResultType>::evaluate(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			arma::Col<ValueType>& result) const
{
	m_fmmTransform.evaluateTrial(point, normal, 
		m_khat, m_nodeCentre, m_nodeSize, result);
}

// Test function

template <typename ResultType>
FmmFunctionMultiplyingTest<ResultType>::FmmFunctionMultiplyingTest(
		const arma::Col<CoordinateType>& khat, 
		const arma::Col<CoordinateType>& nodeCentre,
		const arma::Col<CoordinateType>& nodeSize,
		const FmmTransform<ResultType>& fmmTransform)
	 :	m_khat(khat), m_nodeCentre(nodeCentre), m_nodeSize(nodeSize), 
		m_fmmTransform(fmmTransform)
{
}

template <typename ResultType>
int FmmFunctionMultiplyingTest<ResultType>::argumentDimension() const
{
	return 3;
}

template <typename ResultType>
int FmmFunctionMultiplyingTest<ResultType>::resultDimension() const 
{
	return 1;
}

template <typename ResultType>
inline void FmmFunctionMultiplyingTest<ResultType>::evaluate(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			arma::Col<ValueType>& result) const
{
	m_fmmTransform.evaluateTest(point, normal, 
		m_khat, m_nodeCentre, m_nodeSize, result);
}


// FmmSingleLayerHighFreqFMM

template <typename ValueType>
void FmmSingleLayerHighFreq<ValueType>::evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = nodeCentre - point;
	result(0) =  exp( -kappa*dot(khat, r) );
}

template <typename ValueType>
void FmmSingleLayerHighFreq<ValueType>::evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = point - nodeCentre;
	result(0) =  exp( -kappa*dot(khat, r) );
}

// FmmDoubleLayerHighFreq

template <typename ValueType>
void FmmDoubleLayerHighFreq<ValueType>::evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = nodeCentre - point;
	result(0) =  kappa*exp( -kappa*dot(khat, r) )*dot(khat, normal);
}

template <typename ValueType>
void FmmDoubleLayerHighFreq<ValueType>::evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = point - nodeCentre;
	result(0) =  exp( -kappa*dot(khat, r) );
}

// FmmAdjointDoubleLayerHighFreq

template <typename ValueType>
void FmmAdjointDoubleLayerHighFreq<ValueType>::evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = nodeCentre - point;
	result(0) =  exp( -kappa*dot(khat, r) );
}

template <typename ValueType>
void FmmAdjointDoubleLayerHighFreq<ValueType>::evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = point - nodeCentre;
	result(0) =  -kappa*exp( -kappa*dot(khat, r) )*dot(khat, normal);
}

// FmmHypersingularHighFreq

// The translational error in moving between levels in the octree increases from 
// the single layer potential to the double layer potential to the hypersingular
// operators. This means effectively, only a few octree levels can be used for the 
// hypersingular operator. The situtation is improved a little by using a finer mesh.

// Note that the results are only valid for P1 trial functions. If, for the sake of 
// consistency, we want the ACA and FMM to give the same results for P0 trial functions
// some work must be done. For P0 trial functions, the far field is very similar
// to the P1 solution. The near field where the major difference lies. This same 
// difference is seen with the ACA for P0 trial functions. The hybrid assembly mode
// gives different results to the global and local modes (which are identical).

template <typename ValueType>
void FmmHypersingularHighFreq<ValueType>::evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = nodeCentre - point;
	result(0) =  kappa*exp( -kappa*dot(khat, r) )*dot(khat, normal);
}

template <typename ValueType>
void FmmHypersingularHighFreq<ValueType>::evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = point - nodeCentre;
	result(0) =  kappa*exp( -kappa*dot(khat, r) )*dot(khat, normal);
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmDoubleLayerHighFreq);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmAdjointDoubleLayerHighFreq);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmSingleLayerHighFreq);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmHypersingularHighFreq);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmFunctionMultiplyingTrial);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmFunctionMultiplyingTest);

} // namespace Bempp
