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


#ifndef bempp_grid_function_hpp
#define bempp_grid_function_hpp

#include "../common/common.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/shared_ptr.hpp"

#include "../grid/vtk_writer.hpp"
#include "../fiber/quadrature_strategy.hpp"
#include "../fiber/scalar_traits.hpp"

#include <boost/mpl/set.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/utility/enable_if.hpp>
#include <memory>

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ValueType> class Basis;
template <typename ResultType> class LocalAssemblerForGridFunctions;
template <typename ValueType> class Function;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \cond FORWARD_DECL */
class AssemblyOptions;
class GeometryFactory;
class Grid;
template <int codim> class Entity;
template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType> class Context;
/** \cond endcond */

using Fiber::Function;

/** \ingroup assembly_functions
 *  \brief Function defined on a grid.
 *
 *  This class represents a function defined on a grid and expanded in a
 *  particular function space. */
template <typename BasisFunctionType, typename ResultType>
class GridFunction
{
public:
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;
    typedef typename Fiber::ScalarTraits<ResultType>::RealType MagnitudeType;

    enum DataType { COEFFICIENTS, PROJECTIONS };

    /** \brief Constructor.
     *
     *  Construct an uninitialized grid function. The only way to
     *  initialize it later is using the assignment operator. */
    GridFunction();

    /** \brief Constructor.
     *
     *  \param[in] context      Assembly context from which a quadrature
     *                          strategy can be retrieved.
     *  \param[in] space        Function space to expand the grid function in.
     *  \param[in] dualSpace    Function space dual to \p space.
     *  \param[in] data
     *    If <tt>dataType == COEFFICIENTS</tt>, the vector \p data should have
     *    length <tt>space.globalDofCount()</tt> and contain the expansion
     *    coefficients of the grid function in the space \p space. Otherwise,
     *    if <tt>dataType == PROJECTIONS</tt>, \p data should have length
     *    <tt>dualSpace.globalDofCount()</tt> contain the scalar products of
     *    the grid function and the basis functions of the space \p dualSpace.
     *  \param[in] dataType     Interpretation of the vector
     *                          passed via the argument \p data.
     *
     *  This constructor builds a grid function from either its coefficients in
     *  a function space or from its projections on the basis functions of
     *  another (dual) function space. If the other type of data turns out to
     *  be needed later (for example, the projections vector if coefficients
     *  were supplied in the constructor), it is calculated automatically. The
     *  Context object given in the constructor is then used to determine the
     *  strategy for evaluating any necessary integrals.
     *
     *  \p space and \p dualSpace must be defined on the same grid. */
    GridFunction(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& space,
                 const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
                 const arma::Col<ResultType>& data,
                 DataType dataType);

    /** \brief Constructor.
     *
     *  \param[in] context      Assembly context from which a quadrature
     *                          strategy can be retrieved.
     *  \param[in] space        Function space to expand the grid function in.
     *  \param[in] dualSpace    Function space dual to \p space.
     *  \param[in] coefficients Vector of length <tt>space.globalDofCount()</tt>
     *                          containing the expansion coefficients of the grid
     *                          function in the space \p space.
     *  \param[in] projections  Vector of length <tt>dualSpace.globalDofCount()</tt>
     *                          containing the scalar products of the grid
     *                          function and the basis functions of the space
     *                          \p dualSpace.
     *
     *  \p space and \p dualSpace must be defined on the same grid.
     *
     *  \note This constructor is mainly intended for internal use in BEM++. */
    GridFunction(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& space,
                 const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
                 const arma::Col<ResultType>& coefficients,
                 const arma::Col<ResultType>& projections);

    /** \brief Constructor.
     *
     *  \param[in] context      Assembly context from which a quadrature
     *                          strategy can be retrieved.
     *  \param[in] space        Function space to expand the grid function in.
     *  \param[in] dualSpace    Function space dual to \p space.
     *  \param[in] function     Function object whose values on
     *                          <tt>space.grid()</tt> will be used to construct
     *                          the new grid function.
     *
     *  This constructor builds a grid function by approximating the function
     *  \p function in the basis of space \p space.
     *
     *  \p space and \p dualSpace must be defined on the same grid. */
    GridFunction(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& space,
                 const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
                 const Function<ResultType>& function);

    /** \brief Return whether this function has been properly initialized. */
    bool isInitialized() const;

    /** \brief Grid on which this function is defined.
     *
     * \note An exception is thrown if this function is called on an
     * uninitialized GridFunction object. */
    shared_ptr<const Grid> grid() const;

    /** \brief Space in which this function is expanded. */
    shared_ptr<const Space<BasisFunctionType> > space() const;

    /** \brief Space dual to the space in which this function is expanded. */
    shared_ptr<const Space<BasisFunctionType> > dualSpace() const;

    /** \brief Assembly context used to retrieve the strategy for evaluating
     *  any necessary integrals. */
    shared_ptr<const Context<BasisFunctionType, ResultType> > context() const;

    /** \brief Number of components of this function.
     *
     * \note An exception is thrown if this function is called on an
     * uninitialized GridFunction object. */
    int componentCount() const;

    /** \brief Vector of expansion coefficients of this function in the basis
     *  of its primal space.
     *
     *  If the grid function was constructed from the vector of its projections
     *  on the basis of its dual space, the expansion coefficients in the
     *  primal space are calculated automatically on the first call to
     *  coefficients() and cached internally before being returned.
     *
     *  An exception is thrown if this function is called on an uninitialized
     *  GridFunction object (one constructed with the default constructor). */
    const arma::Col<ResultType>& coefficients() const;

    /** \brief Vector of scalar products of this function with the basis
     *  functions of its dual space.
     *
     *  If the grid function was constructed from the vector of its expansion
     *  coefficients in its primal space, the projections on the basis
     *  functions of its dual space are calculated automatically on the first
     *  call to projections() and cached internally before being returned.
     *
     *  An exception is thrown if this function is called on an uninitialized
     *  GridFunction object (one constructed with the default constructor). */
    const arma::Col<ResultType>& projections() const;

    /** \brief Reset the expansion coefficients of this function in the basis
     *  of its primal space.
     *
     *  As a side effect, any internally stored vector of the projections of
     *  this grid function on the basis functions of its dual space is marked as
     *  invalid and recalculated on the next call to projections(). */
    void setCoefficients(const arma::Col<ResultType>& coeffs);

    /** \brief Reset the vector of scalar products of this function with the basis
     *  functions of its dual space.
     *
     *  As a side effect, any internally stored vector of the coefficients of
     *  this grid function in its primal space is marked as invalid and
     *  recalculated on the next call to coefficients(). */
    void setProjections(const arma::Col<ResultType>& projects);

    /** \brief Return the \f$L^2\f$-norm of the grid function. */
    MagnitudeType L2Norm() const;

    const Fiber::Basis<BasisFunctionType>& basis(const Entity<0>& element) const;

    /** \brief Retrieve the expansion coefficients of this function on a single element.
     *
     *  \param[in] element An element belonging to the grid <tt>space.grid()</tt>.
     *  \param[out] coeffs Vector of the expansion coefficients of this function
     *                     corresponding to the basis functions of the primal space
     *                     living on element \p element.
     *
     *  \note The results of calling this function on an uninitialized
     *  GridFunction object are undefined. */
    void getLocalCoefficients(const Entity<0>& element,
                              std::vector<ResultType>& coeffs) const;

    /** \brief Export this function to a VTK file.

      \param[in] dataType
        Determines whether data are attaches to vertices or cells.

      \param[in] dataLabel
        Label used to identify the function in the VTK file.

      \param[in] fileNamesBase
        Base name of the output files. It should not contain any directory
        part or filename extensions.

      \param[in] filesPath
        Output directory. Can be set to NULL, in which case the files are
        output in the current directory.

      \param[in] type
        Output type (default: ASCII). See Dune reference manual for more
        details.

      \note An exception is thrown if this function is called on an
        uninitialized GridFunction object. */
    void exportToVtk(VtkWriter::DataType dataType,
                     const char* dataLabel,
                     const char* fileNamesBase, const char* filesPath = 0,
                     VtkWriter::OutputType type = VtkWriter::ASCII) const;

    /** \brief Evaluate function at either vertices or barycentres.
     *
     *  \note The results of calling this function on an uninitialized
     *  GridFunction object are undefined. */
    void evaluateAtSpecialPoints(
            VtkWriter::DataType dataType, arma::Mat<ResultType>& result_) const;

private:
    shared_ptr<const Context<BasisFunctionType, ResultType> > m_context;
    shared_ptr<const Space<BasisFunctionType> > m_space;
    shared_ptr<const Space<BasisFunctionType> > m_dualSpace;
    mutable shared_ptr<const arma::Col<ResultType> > m_coefficients;
    mutable shared_ptr<const arma::Col<ResultType> > m_projections;
};

// Overloaded operators

/** \brief Return a copy of the passed function \p g. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator+(
        const GridFunction<BasisFunctionType, ResultType>& g);

/** \brief Return the grid function representing the function \p g
 *  multipled by -1.
 *
 *  \note An exception is thrown if \p g is uninitialized. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator-(
        const GridFunction<BasisFunctionType, ResultType>& g);

/** \brief Return a grid function representing the sum of the operands.
 *
 *  \note Both operands must be initialized and their primal and dual
 *  spaces must be identical, otherwise an exception is thrown. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator+(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const GridFunction<BasisFunctionType, ResultType>& g2);

/** \brief Return a grid function representing the difference of the operands.
 *
 *  \note Both operands must be initialized and their primal and dual
 *  spaces must be identical, otherwise an exception is thrown. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator-(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const GridFunction<BasisFunctionType, ResultType>& g2);

/** \brief Return the grid function representing the function \p g
 *  multiplied by \p scalar.
 *
 *  \note An exception is thrown if \p g is uninitialized. */
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const GridFunction<BasisFunctionType, ResultType>& g, const ScalarType& scalar);

// This type machinery is needed to disambiguate between this operator and
// the one taking a AbstractBoundaryOperator and a GridFunction
/** \brief Return the grid function representing the function \p g
 *  multiplied by \p scalar.
 *
 *  \note An exception is thrown if \p g is uninitialized. */
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    GridFunction<BasisFunctionType, ResultType> >::type
operator*(
        const ScalarType& scalar, const GridFunction<BasisFunctionType, ResultType>& g);

/** \brief Return the grid function representing the function \p g
 *  divided by \p scalar.
 *
 *  \note An exception is thrown if \p g is uninitialized. */
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator/(
        const GridFunction<BasisFunctionType, ResultType>& g1, const ScalarType& scalar);

} // namespace Bempp

#endif
