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

#include "vector.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/lazy.hpp"
#include "../common/shared_ptr.hpp"

#include "../grid/vtk_writer.hpp"
#include "../fiber/quadrature_strategy.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../fiber/surface_normal_dependent_function.hpp"
#include "../fiber/surface_normal_independent_function.hpp"

#include <boost/mpl/set.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/utility/enable_if.hpp>
#include <memory>

namespace Fiber
{

template <typename ValueType> class Basis;
template <typename ResultType> class LocalAssemblerForGridFunctions;
template <typename ValueType> class Function;

} // namespace Fiber

namespace Bempp
{

class AssemblyOptions;
class GeometryFactory;
class Grid;
template <int codim> class Entity;
template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType> class Context;

using Fiber::Function;

/** \brief Function defined on a grid.
 *  \ingroup assembly
 */
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
     *  initialize it later is using the assigment operator. */
    GridFunction();

    // TODO: clarify this description
    /** \brief Constructor.
     *
     * \param[in] space        Function space to expand the grid function in.
     * \param[in] dualSpace    Function space dual to \p space.
     * \param[in] data
     *   Depending on the value of \dataType, vector of the expansion
     *   coefficients of the grid function in the space \p space. or vector of
     *   the scalar products of the grid function and the basis functions of
     *   the space \p dualSpace (in other words, the element of \p dualSpace
     *   yielded by the Riesz map of \p space applied to this grid function).
     * \param[in] dataType     Interpretation of the vector
     *                         passed via the argument \p data.
     *
     * \note End users should not need to call this constructor directly.
     * Use instead one of the "non-member constructors"
     * <tt>gridFunctionFrom...()</tt>. */
    GridFunction(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& space,
                 const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
                 const arma::Col<ResultType>& data,
                 DataType dataType);

    /** \brief Constructor.
     *
     * \param[in] space        Function space to expand the grid function in.
     * \param[in] dualSpace    Function space dual to \p space.
     * \param[in] coefficients Vector of the expansion coefficients of the grid
     *                         function in the space \p space.
     * \param[in] projections  Vector of the scalar products of the grid
     *                         function and the basis functions of the space
     *                         \p dualSpace (in other words, the element of
     *                         \p dualSpace yielded by the Riesz map of
     *                         \p space applied to this grid function).
     *
     * \note End users should not need to call this constructor directly.
     * Use instead one of the "non-member constructors"
     * <tt>gridFunctionFrom...()</tt>. */
    GridFunction(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& space,
                 const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
                 const arma::Col<ResultType>& coefficients,
                 const arma::Col<ResultType>& projections);

    GridFunction(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& space,
                 const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
                 const Function<ResultType>& function);

    /** \brief True if the grid has been properly initialized. */
    bool isInitialized() const;

    /** \brief Grid on which this function is defined. */
    const Grid& grid() const;

    /** \brief Space in which this function is expanded. */
    shared_ptr<const Space<BasisFunctionType> > space() const;

    /** \brief Space dual to the space in which this function is expanded. */
    shared_ptr<const Space<BasisFunctionType> > dualSpace() const;

    shared_ptr<const Context<BasisFunctionType, ResultType> > context() const;

    int codomainDimension() const;

    const arma::Col<ResultType>& coefficients() const;
    const arma::Col<ResultType>& projections() const;
    void setCoefficients(const arma::Col<ResultType>& coeffs);
    void setProjections(const arma::Col<ResultType>& projects);

    /** \brief Return the \$fL^2\f$-norm of the grid function. */
    MagnitudeType L2Norm() const;

    const Fiber::Basis<BasisFunctionType>& basis(const Entity<0>& element) const;
    void getLocalCoefficients(const Entity<0>& element,
                              std::vector<ResultType>& coeffs) const;

    /** \brief Export the function to a VTK file.

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
        details. */
    void exportToVtk(VtkWriter::DataType dataType,
                     const char* dataLabel,
                     const char* fileNamesBase, const char* filesPath = 0,
                     VtkWriter::OutputType type = VtkWriter::ASCII) const;

    /** \brief Evaluate function at either vertices or barycentres. */
    void evaluateAtSpecialPoints(
            VtkWriter::DataType dataType, arma::Mat<ResultType>& result_) const;

private:
    shared_ptr<const Context<BasisFunctionType, ResultType> > m_context;
    shared_ptr<const Space<BasisFunctionType> > m_space;
    shared_ptr<const Space<BasisFunctionType> > m_dualSpace;
    mutable shared_ptr<const arma::Col<ResultType> > m_coefficients;
    mutable shared_ptr<const arma::Col<ResultType> > m_projections;
    // mutable arma::Col<ResultType> m_coefficients;
    // mutable arma::Col<ResultType> m_projections;
};

// Overloaded operators

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator+(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const GridFunction<BasisFunctionType, ResultType>& g2);

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator-(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const GridFunction<BasisFunctionType, ResultType>& g2);

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const GridFunction<BasisFunctionType, ResultType>& g1, const ScalarType& scalar);

// This type machinery is needed to disambiguate between this operator and
// the one taking a AbstractBoundaryOperator and a GridFunction
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    GridFunction<BasisFunctionType, ResultType> >::type
operator*(
        const ScalarType& scalar, const GridFunction<BasisFunctionType, ResultType>& g2);

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator/(
        const GridFunction<BasisFunctionType, ResultType>& g1, const ScalarType& scalar);

} // namespace Bempp

#endif
