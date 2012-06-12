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

#ifndef bempp_potential_hpp
#define bempp_potential_hpp

#include "../common/common.hpp"

#include "../fiber/local_assembler_factory.hpp"
#include "../common/scalar_traits.hpp"

#include "../common/armadillo_fwd.hpp"
#include <memory>

namespace Bempp
{

class EvaluationOptions;
class GeometryFactory;
class Grid;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename ResultType> class InterpolatedFunction;

/** \brief Potential.
 *
 * This class provides the interface for evaluation of a potential \f$k(x)\f$
 * defined by the formula
 *
 * \f[ k(x) = \int_\Gamma K(x, y) \, g(y) \, \mathrm{d}\Gamma, \f]
 *
 * where the integration goes over a surface \f$\Gamma\f$ and the point \f$x\f$
 * does not lie on \f$\Gamma\f$.
 */
template <typename BasisFunctionType_, typename ResultType_>
class Potential
{
public:
    /** \copydoc LinearOperator::BasisFunctionType */
    typedef BasisFunctionType_ BasisFunctionType;
    /** \copydoc LinearOperator::ResultType */
    typedef ResultType_ ResultType;
    /** \copydoc LinearOperator::CoordinateType */
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
    /** \copydoc LinearOperator::LocalAssemblerFactory. */
    typedef Fiber::LocalAssemblerFactory<BasisFunctionType, ResultType, GeometryFactory>
    LocalAssemblerFactory;

    /** \brief Destructor */
    virtual ~Potential() {}

    /** \brief Evaluate the potential of a given moment on a prescribed grid.
     *
     * \param[in] argument
     *   Argument of the potential (\f$\g(y)\f$ in the notation above),
     *   represented by a grid function.
     * \param[in] evaluationGrid
     *   Surface or volume grid at whose vertices the potential will be evaluated.
     * \param[in] assemblerFactory
     *   Factory that will be used to generate appropriate evaluator objects.
     * \param[in] options
     *   Options.
     *
     * \returns The potential represented by a function interpolated on the
     * vertices of \p evaluationGrid.
     *
     * \note This function is not designed to yield accurate values of the
     * potential on the surface containing the charge distribution, i.e.
     * <tt>argument.grid()</tt>. Hence values of the potential at any vertices
     * of \p evaluationGrid that coincide with <tt>argument.grid()</tt> can be
     * badly wrong.
     */
    virtual std::auto_ptr<InterpolatedFunction<ResultType> > evaluateOnGrid(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const Grid& evaluationGrid,
            const Fiber::LocalAssemblerFactory<
            BasisFunctionType, ResultType, GeometryFactory>& assemblerFactory,
            const EvaluationOptions& options) const = 0;

    virtual arma::Mat<ResultType> evaluateAtPoints(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const arma::Mat<CoordinateType>& evaluationPoints,
            const Fiber::LocalAssemblerFactory<
            BasisFunctionType, ResultType, GeometryFactory>& assemblerFactory,
            const EvaluationOptions& options) const = 0;
};

} // namespace Bempp

#endif
