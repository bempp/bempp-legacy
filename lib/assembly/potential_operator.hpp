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

#ifndef bempp_potential_operator_hpp
#define bempp_potential_operator_hpp

#include "../common/common.hpp"

#include "../fiber/quadrature_strategy.hpp"
#include "../common/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"

#include "../common/armadillo_fwd.hpp"
#include <memory>

namespace Bempp
{

/** \cond FORWARD_DECL */
class EvaluationOptions;
class EvaluationOptions;
class GeometryFactory;
class Grid;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename ResultType> class InterpolatedFunction;
template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType>
class AssembledPotentialOperator;

/** \endcond */

/** \ingroup potential_operators
 *  \brief Potential operator.
 *
 *  This class represents a linear operator that, acting on a function \f$g\f$
 *  defined on a surface \f$\Gamma\f$ embedded in a space \f$\Omega\f$ of
 *  dimension higher by one, produces a *potential* defined at any point of
 *  \f$\Omega\f$ lying outside \f$\Gamma\f$. The function \f$g\f$ is called the
 *  *charge distribution*.
 *
 *  The functions evaluateOnGrid() and evaluateAtPoints() can be used to
 *  evaluate the potential produced by a given charge distribution, represented
 *  with a GridFunction object, at specified points in \f$\Omega \setminus \Gamma\f$.
 *
 *  \tparam BasisFunctionType_
 *    Type of the values of the (components of the) basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam ResultType_
 *    Type of the values of the (components of the) potential.
 *
 *  Both template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 *  Both types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. If \p BasisFunctionType_ is
 *  set to a complex type, then \p ResultType_ must be set to the same type.
 */
template <typename BasisFunctionType_, typename ResultType_>
class PotentialOperator
{
public:
    /** \brief Type of the values of the (components of the) basis functions
     * into which functions acted upon by the operator are expanded. */
    typedef BasisFunctionType_ BasisFunctionType;
    /** \brief Type of the values of the (components of the) potential. */
    typedef ResultType_ ResultType;
    /** \brief Type used to represent coordinates. */
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
    /** \brief Type of the appropriate instantiation of
     *  Fiber::QuadratureStrategy. */
    typedef Fiber::QuadratureStrategy<BasisFunctionType, ResultType, GeometryFactory>
    QuadratureStrategy;
    /** \brief Type of the appropriate instantiation of
     *  Fiber::EvaluatorForIntegralOperators. */
    typedef Fiber::EvaluatorForIntegralOperators<ResultType> Evaluator;

    /** \brief Destructor */
    virtual ~PotentialOperator() {}

    /** \brief Evaluate the potential of a given charge distribution on a
     *  prescribed grid.
     *
     * \param[in] argument
     *   Argument of the potential operator (\f$\psi(y)\f$ in the notation
     *   above), represented by a grid function.
     * \param[in] evaluationGrid
     *   Grid at whose vertices the potential will be evaluated. The grid may
     *   have arbitrary dimension, but must be embedded in a world of the same
     *   dimension as <tt>argument.grid()</tt>.
     * \param[in] quadStrategy
     *   A #QuadratureStrategy object controlling how the integrals will be
     *   evaluated.
     * \param[in] options
     *   Evaluation options.
     *
     * \returns The potential represented by a function interpolated on the
     * vertices of \p evaluationGrid.
     *
     * \note This function is not designed to yield accurate values of the
     * potential on the surface \f$\Gamma\f$ containing the charge
     * distribution, i.e. <tt>argument.grid()</tt>, even if the potential has a
     * unique extension from \f$\Omega \setminus \Gamma\f$ to \f$\Gamma\f$.
     * Hence values of the potential at any vertices of \p evaluationGrid that
     * coincide with \f$\Gamma\f$ can be badly wrong.
     *
     * The current implementation does not yet take special measures to prevent
     * loss of accuracy *near* \f$\Gamma\f$, either. If in doubt, increase the
     * quadrature accuracy. */
    virtual std::auto_ptr<InterpolatedFunction<ResultType> > evaluateOnGrid(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const Grid& evaluationGrid,
            const QuadratureStrategy& quadStrategy,
            const EvaluationOptions& options) const = 0;

    /** \brief Evaluate the potential of a given charge distribution at
     *  prescribed points.
     *
     * \param[in] argument
     *   Argument of the potential operator (\f$\psi(y)\f$ in the notation above),
     *   represented by a grid function.
     * \param[in] evaluationPoints
     *   2D array whose (i, j)th element is the ith coordinate of the jth point
     *   at which the potential should be evaluated. The first dimension of this
     *   array should be equal to <tt>argument.grid().dimWorld()</tt>.
     * \param[in] quadStrategy
     *   A #QuadratureStrategy object controlling how the integrals will be
     *   evaluated.
     * \param[in] options
     *   Evaluation options.
     *
     * \returns A 2D array whose (i, j)th element is the ith component of the
     * potential at the jth point.
     *
     * \note This function is not designed to yield accurate values of the
     * potential on the surface \f$\Gamma\f$ containing the charge
     * distribution, i.e. <tt>argument.grid()</tt>, even if the potential has a
     * unique extension from \f$\Omega \setminus \Gamma\f$ to \f$\Gamma\f$.
     * Hence values of the potential at any points belonging to \f$\Gamma\f$
     * can be badly wrong.
     *
     * The current implementation does not yet take special measures to prevent
     * loss of accuracy *near* \f$\Gamma\f$, either. Users are advised to
     * increase the quadrature accuracy for points lying in the vicinity of
     * \f$\Gamma\f$. */
    virtual arma::Mat<ResultType> evaluateAtPoints(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const arma::Mat<CoordinateType>& evaluationPoints,
            const QuadratureStrategy& quadStrategy,
            const EvaluationOptions& options) const = 0;

    /** \brief Create and return an AssembledPotentialOperator object.
     *
     *  The returned AssembledPotentialOperator object stores the values of the
     *  potentials generated at the points listed in the array \p
     *  evaluationPoints by the charge distributions equal to the individual
     *  basis functions of the space \p space. The object can afterwards be used
     *  to evaluate efficiently the potentials generated by multiple
     *  GridFunctions expanded in the space \p space.
     *
     * \param[in] space
     *   The space whose basis functions will be taken as the charge distributions
     *   inducing the potentials to be evaluated.
     * \param[in] evaluationPoints
     *   2D array whose (i, j)th element is the ith coordinate of the jth point
     *   at which the potential should be evaluated. The first dimension of this
     *   array should be equal to <tt>space.grid().dimWorld()</tt>.
     * \param[in] quadStrategy
     *   A #QuadratureStrategy object controlling how the integrals will be
     *   evaluated.
     * \param[in] options
     *   Evaluation options. This parameter controls, notably, the format
     *   used to store the matrix of precalculated potential values: a dense
     *   matrix or an H-matrix.
     */
    virtual AssembledPotentialOperator<BasisFunctionType, ResultType>
    assemble(
            const shared_ptr<const Space<BasisFunctionType> >& space,
            const shared_ptr<const arma::Mat<CoordinateType> >& evaluationPoints,
            const QuadratureStrategy& quadStrategy,
            const EvaluationOptions& options) const = 0;

    /** \brief Create and return an evaluator object that can be used to optimise
     *  repeated application of this potential operator to a particular grid
     *  function at different points.
     *
     *  Use of the returned evaluator object is advantageous compared to calling
     *  evaluateAtPoints() or evaluateOnGrid() when the evaluation points are
     *  determined iteratively rather than being all known in advance.
     *
     *  Example:
        \code
        GridFunction<double, double> function = ...;
        Laplace3dSingleLayerPotentialOperator<double, double> slPotOp;
        typedef Fiber::EvaluatorForIntegralOperators<double> Evaluator;
        NumericalQuadratureStrategy<double, double> quadStrategy;
        EvaluationOptions evaluationOptions;
        std::auto_ptr<Evaluator> slPotOpEvaluator =
            slPotOp.makeEvaluator(function, quadStrategy, evaluationOptions);
        // Evaluate the potential at points (1, 0, 0) and (2, 3, 4)
        arma::Mat<double> points(3, 2);
        points(0, 0) = 1.; points(1, 0) = 0.; points(2, 0) = 0.;
        points(0, 1) = 2.; points(1, 1) = 3.; points(2, 1) = 4.;
        arma::Mat<double> slPotValues;
        slPotEvaluator->evaluate(Evaluator::FAR_FIELD,
                                 points, slPotValues);
        \endcode
     *
     * \param[in] argument
     *   The grid function to which this potential operator will be applied.
     * \param[in] quadStrategy
     *   A #QuadratureStrategy object controlling how the integrals will be
     *   evaluated.
     * \param[in] options
     *   Evaluation options.
     */
    virtual std::auto_ptr<Evaluator> makeEvaluator(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const QuadratureStrategy& quadStrategy,
            const EvaluationOptions& options) const = 0;

    /** \brief Number of components of the values of the potential.
     *
     *  E.g. 1 for a scalar-valued potential, 3 for a vector-valued potential.
     */
    virtual int componentCount() const = 0;
};

} // namespace Bempp

#endif
