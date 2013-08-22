
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
		const arma::Col<CoordinateType>& centre,
		const FmmTransform<ResultType>& fmmTransform)
	 : m_khat(khat), m_centre(centre), m_fmmTransform(fmmTransform)
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
		m_khat, m_centre, result);
}

// Test function

template <typename ResultType>
FmmFunctionMultiplyingTest<ResultType>::FmmFunctionMultiplyingTest(
		const arma::Col<CoordinateType>& khat, 
		const arma::Col<CoordinateType>& centre,
		const FmmTransform<ResultType>& fmmTransform)
	 : m_khat(khat), m_centre(centre), m_fmmTransform(fmmTransform)
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
		m_khat, m_centre, result);
}

// FmmDoubleLayerHighFreq

template <typename ValueType>
void FmmDoubleLayerHighFreq<ValueType>::evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& centre,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = centre - point;
	result(0) =  kappa*exp( -kappa*dot(khat, r) )*dot(khat, normal);
}

template <typename ValueType>
void FmmDoubleLayerHighFreq<ValueType>::evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& centre,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = point - centre;
	result(0) =  exp( -kappa*dot(khat, r) );
}

// FmmAdjointDoubleLayerHighFreq

template <typename ValueType>
void FmmAdjointDoubleLayerHighFreq<ValueType>::evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& centre,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = centre - point;
	result(0) =  exp( -kappa*dot(khat, r) );
}

template <typename ValueType>
void FmmAdjointDoubleLayerHighFreq<ValueType>::evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& centre,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = point - centre;
	result(0) =  -kappa*exp( -kappa*dot(khat, r) )*dot(khat, normal);
}

// FmmSingleLayerHighFreqFMM

template <typename ValueType>
void FmmSingleLayerHighFreq<ValueType>::evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& centre,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = centre - point;
	result(0) =  exp( -kappa*dot(khat, r) );
}

template <typename ValueType>
void FmmSingleLayerHighFreq<ValueType>::evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& centre,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = point - centre;
	result(0) =  exp( -kappa*dot(khat, r) );
}

// FmmHypersingularHighFreq

template <typename ValueType>
void FmmHypersingularHighFreq<ValueType>::evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& centre,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = centre - point;
	throw std::invalid_argument("FmmHypersingularHighFreq::evaluateTrial(): "
			"FMM not currently implemented for the hypersingular operator");
//	result(0) =  kappa*exp( -kappa*dot(khat, r) )*dot(khat, normal);
}

template <typename ValueType>
void FmmHypersingularHighFreq<ValueType>::evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& khat,
			const arma::Col<CoordinateType>& centre,
			arma::Col<ValueType>& result) const
{
	ValueType kappa = FmmHighFreq<ValueType>::kappa();
	arma::Col<CoordinateType> r = point - centre;
	throw std::invalid_argument("FmmHypersingularHighFreq::evaluateTest(): "
			"FMM not currently implemented for the hypersingular operator");
//	result(0) =  -kappa*exp( -kappa*dot(khat, r) )*dot(khat, normal);
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmDoubleLayerHighFreq);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmAdjointDoubleLayerHighFreq);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmSingleLayerHighFreq);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmHypersingularHighFreq);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmFunctionMultiplyingTrial);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmFunctionMultiplyingTest);

} // namespace Bempp
