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

#ifndef bempp_fmm_integrator_hpp
#define bempp_fmm_integrator_hpp

#include "../common/common.hpp"

#include "../common/types.hpp" // Point3D
#include "../common/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"

#include "../fiber/numerical_quadrature.hpp"
#include "../fiber/test_function_integrator.hpp"

#include <tbb/concurrent_unordered_map.h>
#include <map>
#include <vector>

namespace Fiber
{
/** \cond FORWARD_DECL */
template <typename CoordinateType> class CollectionOfBasisTransformations;
template <typename UserFunctionType> class Function;
class OpenClHandler;
template <typename CoordinateType> class RawGridGeometry;
/** \endcond */
}

namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class Space;
class AssemblyOptions;
class GeometryFactory;
/** \endcond */

template <typename BasisFunctionType, typename ResultType>
class FmmLocalAssembler
{
public:
	typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
	typedef ResultType UserFunctionType; // get operator mismatches if not the same as ResultType

	FmmLocalAssembler(const Space<BasisFunctionType>& space,
		const AssemblyOptions& options, bool conjugateBasis = true);

	void setFunction(Fiber::Function<ResultType> *function);

	void evaluateLocalWeakForms(
            const std::vector<int>& elementIndices,
            std::vector<arma::Col<ResultType> > & result);

	~FmmLocalAssembler();
private:
	typedef Fiber::TestFunctionIntegrator<BasisFunctionType, ResultType> Integrator;
	typedef tbb::concurrent_unordered_map<Fiber::SingleQuadratureDescriptor,
		Integrator*> IntegratorMap;

	void clearIntegratorMap();

	const Integrator& selectIntegrator(int elementIndex);

	const Integrator& getIntegrator(const Fiber::SingleQuadratureDescriptor& index);


	typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
	typedef std::vector<const Fiber::Shapeset<BasisFunctionType>*> ShapesetPtrVector;

	shared_ptr<Fiber::RawGridGeometry<CoordinateType> > m_rawGeometry;
	shared_ptr<GeometryFactory> m_geometryFactory;
//	shared_ptr<Fiber::Function<ResultType> > m_function;
	Fiber::Function<ResultType> *m_function;
	shared_ptr<Fiber::OpenClHandler> m_openClHandler;
	shared_ptr<ShapesetPtrVector> m_shapesets;
	shared_ptr<const Fiber::CollectionOfBasisTransformations<CoordinateType> > 
		m_transformations;

	IntegratorMap m_testFunctionIntegrators;
	bool m_conjugateBasis;
};

} // namespace Bempp

#endif
