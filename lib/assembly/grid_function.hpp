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
#include "../common/deprecated.hpp"
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
/** \endcond */

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

    // Recommended constructors

    /** \brief Constructor.
     *
     *  Construct an uninitialized grid function. The only way to
     *  initialize it later is using the assignment operator. */
    GridFunction();

    /** Constructor.
     *
     *  \param[in] context      Assembly context from which a quadrature
     *                          strategy can be retrieved.
     *  \param[in] space        Function space to expand the grid function in.
     *  \param[in] coefficients
     *    %Vector of length <tt>space.globalDofCount()</tt> containing the expansion
     *    coefficients of the grid function in the space \p space.
     *
     *  This constructor builds a grid function from its coefficients in a
     *  function space.
     */
    GridFunction(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& space,
                 const arma::Col<ResultType>& coefficients);

    /** Constructor.
     *
     *  \param[in] context      Assembly context from which a quadrature
     *                          strategy can be retrieved.
     *  \param[in] space        Function space to expand the grid function in.
     *  \param[in] dualSpace    Function space dual to \p space.
     *  \param[in] projections
     *    %Vector of length <tt>dualSpace.globalDofCount()</tt> containing the
     *    scalar products of the grid function and the basis functions of the
     *    space \p dualSpace.
     *
     *  This constructor builds a grid function expanded in the basis of the
     *  space \p space from the projections of this function on the basis
     *  functions of another space \p dualSpace. Both spaces must be defined on
     *  the same grid. */
    GridFunction(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& space,
                 const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
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
     *  This constructor builds a grid function belonging to the space
     *  \p space and approximating the function \f$f\f$ defined by the
     *  object \p function. The approximate coefficients \f$\{f_i\}\f$
     *  of \p function in the basis \f$\{\phi_i\}\f$ of \p space are
     *  determined by solving the equation
     *  \f[ \sum_j \langle \psi_i, \phi_j \rangle f_i =
     *      \langle \psi_i, f \rangle \f]
     *  in the least-squares sense, with \f$\{\psi_i\}\f$ denoting the
     *  set of basis functions of \p dualSpace. The two spaces must be
     *  defined on the same grid. */
    GridFunction(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& space,
                 const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
                 const Function<ResultType>& function);

    // Deprecated constructors

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
     *  \p space and \p dualSpace must be defined on the same grid.
     *
     *  \deprecated This constructor is deprecated. In new code, use one of
     *  the other constructors. */
    BEMPP_DEPRECATED
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
     *  \param[in] coefficients %Vector of length <tt>space.globalDofCount()</tt>
     *                          containing the expansion coefficients of the grid
     *                          function in the space \p space.
     *  \param[in] projections  %Vector of length <tt>dualSpace.globalDofCount()</tt>
     *                          containing the scalar products of the grid
     *                          function and the basis functions of the space
     *                          \p dualSpace.
     *
     *  \p space and \p dualSpace must be defined on the same grid.
     *
     *  \deprecated This constructor is deprecated. In new code, use one of
     *  the other constructors. */
    BEMPP_DEPRECATED
    GridFunction(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& space,
                 const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
                 const arma::Col<ResultType>& coefficients,
                 const arma::Col<ResultType>& projections);

    // Member functions

    /** \brief Return whether this function has been properly initialized. */
    bool isInitialized() const;

    /** \brief Grid on which this function is defined.
     *
     * \note An exception is thrown if this function is called on an
     * uninitialized GridFunction object. */
    shared_ptr<const Grid> grid() const;

    /** \brief Space in which this function is expanded. */
    shared_ptr<const Space<BasisFunctionType> > space() const;

    /** \brief Space dual to the space in which this function is expanded.
     *
     *  \deprecated This function is provided only for backward compatibility
     *  reasons. It returns the shared pointer to the dual space of the
     *  GridFunction, if one was supplied during the construction of the
     *  latter, or a null pointer otherwise. */
    shared_ptr<const Space<BasisFunctionType> > dualSpace() const;

    /** \brief Assembly context used to retrieve the strategy for evaluating
     *  any necessary integrals. */
    shared_ptr<const Context<BasisFunctionType, ResultType> > context() const;

    /** \brief Number of components of this function.
     *
     * \note An exception is thrown if this function is called on an
     * uninitialized GridFunction object. */
    int componentCount() const;

    /** \brief %Vector of expansion coefficients of this function in the basis
     *  of its expansion space.
     *
     *  An exception is thrown if this function is called on an uninitialized
     *  GridFunction object (one constructed with the default constructor). */
    const arma::Col<ResultType>& coefficients() const;

    /** \brief %Vector of scalar products of this function with the basis
     *  functions of \p dualSpace.
     *
     *  \p dualSpace must be defined on the same grid as the space in which the
     *  GridFunction is expanded.
     *
     *  An exception is thrown if this function is called on an uninitialized
     *  GridFunction object (one constructed with the default constructor). */
    arma::Col<ResultType> projections(
            const Space<BasisFunctionType>& dualSpace_) const;

    /** \brief %Vector of scalar products of this function with the basis
     *  functions of its dual space.
     *
     *  \deprecated This function is provided only for backward compatibility
     *  purposes. It works only if the dual space was specified during the
     *  construction of the GridFunction. In new code the other overload of
     *  projections() should be used.
     *
     *  An exception is thrown if this function is called on an uninitialized
     *  GridFunction object (one constructed with the default constructor). */
    BEMPP_DEPRECATED arma::Col<ResultType> projections() const;

    /** \brief Reset the expansion coefficients of this function in the basis
     *  of its primal space.
     *
     *  As a side effect, any internally stored vector of the projections of
     *  this grid function on the basis functions of its dual space is marked as
     *  invalid and recalculated on the next call to projections(). */
    void setCoefficients(const arma::Col<ResultType>& coeffs);

    /** \brief Reinitialize the function by specifying the vector of its scalar
     *  products with the basis functions of \p dualSpace.
     *
     *  \p dualSpace must be defined on the same grid as the space in which the
     *  GridFunction is expanded. */
    void setProjections(
            const Space<BasisFunctionType>& dualSpace_,
            const arma::Col<ResultType>& projects);

    /** \brief Reset the vector of scalar products of this function with the
     *  basis functions of its dual space.
     *
     *  \deprecated This function is provided only for backward compatibility
     *  purposes. It works only if the dual space was specified during the
     *  construction of the GridFunction. In new code the other overload of
     *  setProjections() should be used. */
    BEMPP_DEPRECATED void setProjections(const arma::Col<ResultType>& projects);

    /** \brief Return the \f$L^2\f$-norm of the grid function.
     *
     *  \note For better accuracy, prefer to use L2NormOfDifference() or
     *  estimateL2Error() to calculate the norm of the difference between a grid
     *  function and a function defined by an analytical expression.
     */
    MagnitudeType L2Norm() const;

    const Fiber::Basis<BasisFunctionType>& basis(const Entity<0>& element) const;

    /** \brief Retrieve the expansion coefficients of this function on a single element.
     *
     *  \param[in] element An element belonging to the grid <tt>space.grid()</tt>.
     *  \param[out] coeffs %Vector of the expansion coefficients of this function
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
        uninitialized GridFunction object.

      \deprecated This member function is deprecated. Use the standalone
      function ::exportToVtk() instead.
      */
    BEMPP_DEPRECATED
    void exportToVtk(VtkWriter::DataType dataType,
                     const char* dataLabel,
                     const char* fileNamesBase, const char* filesPath = 0,
                     VtkWriter::OutputType type = VtkWriter::ASCII) const;

    /** \brief Export this function to a Gmsh (.msh) file.

      \param[in] dataLabel
        Label used to identify the function in the Gmsh file.

      \param[in] fileName
        Name of the file to be written.

      \note An exception is thrown if this function is called on an
        uninitialized GridFunction object.

      \note Currently this function treats all elements as independent from each
      other, so that for example each internal mesh node is output as many
      times as there are elements adjacent to it. This makes it possible to
      export GridFunctions expanded in spaces with discontinuous bases. However,
      it also means that a \c .msh file generated by exportToGmsh cannot be
      imported back by BEM++ and used to create a Grid object.

      \deprecated This member function is deprecated. Use the standalone
      function ::exportToGmsh() instead.
    */
    BEMPP_DEPRECATED
    void exportToGmsh(
        const char* dataLabel, const char* fileName) const;

    /** \brief Evaluate function at either vertices or barycentres.
     *
     *  \note The results of calling this function on an uninitialized
     *  GridFunction object are undefined. */
    void evaluateAtSpecialPoints(
            VtkWriter::DataType dataType, arma::Mat<ResultType>& result_) const;
    void evaluateAtSpecialPoints(
            VtkWriter::DataType dataType,
            arma::Mat<CoordinateType>& points, arma::Mat<ResultType>& values) const;

   /** \brief Evaluate function at specific points lying on a given element.
    *
    *  \param[in] element   An element belonging to the grid on which \p function
    *                       is defined.
    *  \param[in] local     A 2D array whose (i,j)th element is the ith coordinate of
    *                       the jth point at which the grid function should be
    *                       evaluated, in the local coordinate system of
    *                       \p element.
    *  \param[in] values    A 2D array whose (i,j)th element is the ith component
    *                       of \p function evaluated at the jth point.
    */
    void evaluate(
        const Entity<0>& element,
        const arma::Mat<CoordinateType>& local,
        arma::Mat<ResultType>& values) const;

private:
    void initializeFromCoefficients(
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& space,
            const arma::Col<ResultType>& coefficients);
    void initializeFromProjections(
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const shared_ptr<const Space<BasisFunctionType> >& space,
            const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
            const arma::Col<ResultType>& projections);

private:
    shared_ptr<const Context<BasisFunctionType, ResultType> > m_context;
    shared_ptr<const Space<BasisFunctionType> > m_space;
    shared_ptr<const Space<BasisFunctionType> > m_dualSpace;
    mutable shared_ptr<const arma::Col<ResultType> > m_coefficients;
    // mutable shared_ptr<const arma::Col<ResultType> > m_projections;
};

// Overloaded operators

/** \relates GridFunction
 *  \brief Return a copy of the passed function \p g. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator+(
        const GridFunction<BasisFunctionType, ResultType>& g);

/** \relates GridFunction
 *  \brief Return the grid function representing the function \p g
 *  multipled by -1.
 *
 *  \note An exception is thrown if \p g is uninitialized. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator-(
        const GridFunction<BasisFunctionType, ResultType>& g);

/** \relates GridFunction
 *  \brief Return a grid function representing the sum of the operands.
 *
 *  \note Both operands must be initialized and their primal and dual
 *  spaces must be identical, otherwise an exception is thrown. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator+(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const GridFunction<BasisFunctionType, ResultType>& g2);

/** \relates GridFunction
 *  \brief Return a grid function representing the difference of the operands.
 *
 *  \note Both operands must be initialized and their primal and dual
 *  spaces must be identical, otherwise an exception is thrown. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator-(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const GridFunction<BasisFunctionType, ResultType>& g2);

/** \relates GridFunction
 *  \brief Return the grid function representing the function \p g
 *  multiplied by \p scalar.
 *
 *  \note An exception is thrown if \p g is uninitialized. */
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const GridFunction<BasisFunctionType, ResultType>& g, const ScalarType& scalar);

// This type machinery is needed to disambiguate between this operator and
// the one taking a AbstractBoundaryOperator and a GridFunction
/** \relates GridFunction
 *  \brief Return the grid function representing the function \p g
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

/** \relates GridFunction
 *  \brief Return the grid function representing the function \p g
 *  divided by \p scalar.
 *
 *  \note An exception is thrown if \p g is uninitialized. */
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator/(
        const GridFunction<BasisFunctionType, ResultType>& g1, const ScalarType& scalar);

// Export

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
template <typename BasisFunctionType, typename ResultType>
void exportToVtk(const GridFunction<BasisFunctionType, ResultType>& gridFunction,
                 VtkWriter::DataType dataType,
                 const char* dataLabel,
                 const char* fileNamesBase, const char* filesPath = 0,
                 VtkWriter::OutputType type = VtkWriter::ASCII);

/** \brief Export this function to a Gmsh (.msh) file.

  \param[in] dataLabel
    Label used to identify the function in the Gmsh file.

  \param[in] fileName
    Name of the file to be written.

  \note An exception is thrown if this function is called on an
    uninitialized GridFunction object.

  \note Currently this function treats all elements as independent from each
  other, so that for example each internal mesh node is output as many
  times as there are elements adjacent to it. This makes it possible to
  export GridFunctions expanded in spaces with discontinuous bases. However,
  it also means that a \c .msh file generated by exportToGmsh cannot be
  imported back by BEM++ and used to create a Grid object.
*/
template <typename BasisFunctionType, typename ResultType>
void exportToGmsh(const GridFunction<BasisFunctionType, ResultType>& gridFunction,
                  const char* dataLabel, const char* fileName);

} // namespace Bempp

#endif
