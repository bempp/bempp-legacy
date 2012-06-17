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

#include "../grid/vtk_writer.hpp"
#include "../fiber/local_assembler_factory.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../fiber/surface_normal_dependent_function.hpp"
#include "../fiber/surface_normal_independent_function.hpp"
#include "vector.hpp"

#include "surface_normal_dependent_functor.hpp"
#include "surface_normal_independent_functor.hpp"

#include "../common/armadillo_fwd.hpp"
#include <boost/mpl/set.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/utility/enable_if.hpp>
#include <memory>

namespace Fiber
{

template <typename ValueType> class Basis;
template <typename ValueType> class Function;
template <typename ResultType> class LocalAssemblerForGridFunctions;

} // namespace Fiber

namespace Bempp
{

class AssemblyOptions;
class GeometryFactory;
class Grid;
template <int codim> class Entity;
template <typename BasisFunctionType> class Space;


/** \brief Function defined on a grid.
 *  \ingroup assembly
 */
template <typename BasisFunctionType, typename ResultType>
class GridFunction
{
public:
    typedef Fiber::LocalAssemblerFactory<BasisFunctionType, ResultType, GeometryFactory>
    LocalAssemblerFactory;
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

    /** \brief Constructor.
     *
     * \param[in] space        Function space to expand the grid function in.
     * \param[in] coefficients Vector of the expansion coefficients of the grid
     *                         function in the chosen space.
     * \param[in] projections  Vector of the scalar products of the grid
     *                         function and the basis functions from the chosen
     *                         space.

     * \note End users should not need to call this constructor directly.
     * Use instead one of the "non-member constructors"
     * gridFunctionFrom...(). */
    GridFunction(const Space<BasisFunctionType>& space,
                 const arma::Col<ResultType>& coefficients,
                 const arma::Col<ResultType>& projections);

    /** \brief Grid on which this function is defined. */
    const Grid& grid() const;

    /** \brief Space in which this function is expanded. */
    const Space<BasisFunctionType>& space() const;

    int codomainDimension() const;

    const arma::Col<ResultType>& coefficients() const;
    const arma::Col<ResultType>& projections() const;
    void setCoefficients(const arma::Col<ResultType>& coeffs);
    void setProjections(const arma::Col<ResultType>& projects);

    const Fiber::Basis<BasisFunctionType>& basis(const Entity<0>& element) const;
    void getLocalCoefficients(const Entity<0>& element,
                              std::vector<ResultType>& coeffs) const;

    /** Export the function to a VTK file.

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

private:
    /** \brief Evaluate function at either vertices or barycentres. */
    void evaluateAtSpecialPoints(
            VtkWriter::DataType dataType, arma::Mat<ResultType>& result) const;

private:
    const Space<BasisFunctionType>& m_space;
    arma::Col<ResultType> m_coefficients;
    arma::Col<ResultType> m_projections;
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




template <typename BasisFunctionType>
GridFunction<BasisFunctionType,float>
operator*(const float scalar, const GridFunction<BasisFunctionType, float>& g2)
{
  return g2 * scalar;
}

template <typename BasisFunctionType>
GridFunction<BasisFunctionType,float>
operator*(const double scalar, const GridFunction<BasisFunctionType, float>& g2)
{
  return g2 * scalar;
}


template <typename BasisFunctionType>
GridFunction<BasisFunctionType,double>
operator*(const float scalar, const GridFunction<BasisFunctionType, double>& g2)
{
  return g2 * scalar;
}

template <typename BasisFunctionType>
GridFunction<BasisFunctionType,double>
operator*(const double scalar, const GridFunction<BasisFunctionType, double>& g2)
{
  return g2 * scalar;
}

template <typename BasisFunctionType>
GridFunction<BasisFunctionType,std::complex<float> >
operator*(const float scalar, const GridFunction<BasisFunctionType, std::complex<float> >& g2)
{
  return g2 * scalar;
}

template <typename BasisFunctionType>
GridFunction<BasisFunctionType,std::complex<float> >
operator*(const double scalar, const GridFunction<BasisFunctionType, std::complex<float> >& g2)
{
  return g2 * scalar;
}

template <typename BasisFunctionType>
GridFunction<BasisFunctionType,std::complex<float> >
operator*(const std::complex<float> scalar, const GridFunction<BasisFunctionType, std::complex<float> >& g2)
{
  return g2 * scalar;
}

template <typename BasisFunctionType>
GridFunction<BasisFunctionType,std::complex<float> >
operator*(const std::complex<double> scalar, const GridFunction<BasisFunctionType, std::complex<float> >& g2)
{
  return g2 * scalar;
}

template <typename BasisFunctionType>
GridFunction<BasisFunctionType,std::complex<double> >
operator*(const float scalar, const GridFunction<BasisFunctionType, std::complex<double> >& g2)
{
  return g2 * scalar;
}

template <typename BasisFunctionType>
GridFunction<BasisFunctionType,std::complex<double> >
operator*(const double scalar, const GridFunction<BasisFunctionType, std::complex<double> >& g2)
{
  return g2 * scalar;
}

template <typename BasisFunctionType>
GridFunction<BasisFunctionType,std::complex<double> >
operator*(const std::complex<float> scalar, const GridFunction<BasisFunctionType, std::complex<double> >& g2)
{
  return g2 * scalar;
}

template <typename BasisFunctionType>
GridFunction<BasisFunctionType,std::complex<double> >
operator*(const std::complex<double> scalar, const GridFunction<BasisFunctionType, std::complex<double> >& g2)
{
  return g2 * scalar;
}


// This type machinery is needed to disambiguate between this operator and
// the one taking a LinearOperator and a GridFunction
//template <typename BasisFunctionType, typename ResultType, typename ScalarType>
//typename boost::enable_if<
//    typename boost::mpl::has_key<
//        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
//        ScalarType>,
//    GridFunction<BasisFunctionType, ResultType> >::type
//operator*(
//        const ScalarType& scalar, const GridFunction<BasisFunctionType, ResultType>& g2);

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator/(
        const GridFunction<BasisFunctionType, ResultType>& g1, const ScalarType& scalar);

// Non-member constructors

/** \brief Construct a grid function from its expansion coefficients in a function space.
 *
 * \param[in] space        Function space to expand the grid function in.
 * \param[in] coefficients Expansion coefficients of the grid
 *                         function in the chosen space. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
gridFunctionFromCoefficients(const Space<BasisFunctionType>& space,
                             const arma::Col<ResultType>& coefficients);

/** \brief Construct a grid function from its expansion coefficients in a function space.
 *
 * \param[in] space        Function space to expand the grid function in.
 * \param[in] projections  Scalar products of the grid function and the basis
 *                         functions of the chosen space. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
gridFunctionFromProjections(const Space<BasisFunctionType>& space,
                            const arma::Col<ResultType>& projections);

/** \brief Construct a grid function from a Fiber::Function object. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
gridFunctionFromFiberFunction(
        const Space<BasisFunctionType>& space,
        const Fiber::Function<ResultType>& function,
        const typename GridFunction<BasisFunctionType, ResultType>::LocalAssemblerFactory& factory,
        const AssemblyOptions& assemblyOptions);

/** \brief Construct a grid function from a functor independent from surface orientation.
 *
 * \param[in] space           Function space to expand the grid function in.
 * \param[in] functor         Functor object fulfilling the requirements
 *                            described in the documentation of
 *                            Fiber::SurfaceNormalIndependentFunction.
 * \param[in] factory         Factory to be used to generate the necessary local
 *                            assembler, which will be employed to evaluate
 *                            the scalar products of the grid function with the
 *                            basis functions of the chosen space.
 * \param[in] assemblyOptions Options controlling the assembly procedure.
 *
 * \returns The constructed grid function. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
gridFunctionFromSurfaceNormalIndependentFunctor(
        const Space<BasisFunctionType>& space,
        const SurfaceNormalIndependentFunctor<ResultType>& functor,
        const typename GridFunction<BasisFunctionType, ResultType>::LocalAssemblerFactory& factory,
        const AssemblyOptions& assemblyOptions)
{
    Fiber::SurfaceNormalIndependentFunction<SurfaceNormalIndependentFunctor<ResultType> > fiberFunction(functor);
    return gridFunctionFromFiberFunction(
                space, fiberFunction, factory, assemblyOptions);
}

/** \brief Construct a grid function from a functor dependent on surface orientation.
 *
 * \param[in] space           Function space to expand the grid function in.
 * \param[in] functor         Functor object fulfilling the requirements
 *                            described in the documentation of
 *                            Fiber::SurfaceNormalDependentFunction.
 * \param[in] factory         Factory to be used to generate the necessary local
 *                            assembler, which will be employed to evaluate
 *                            the scalar products of the grid function with the
 *                            basis functions of the chosen space.
 * \param[in] assemblyOptions Options controlling the assembly procedure.
 *
 * \returns The constructed grid function. */
template <typename BasisFunctionType, typename Functor>
GridFunction<BasisFunctionType, typename Functor::ValueType>
gridFunctionFromSurfaceNormalDependentFunctor(
        const Space<BasisFunctionType>& space,
        const Functor& functor,
        const typename GridFunction<BasisFunctionType, typename Functor::ValueType>::LocalAssemblerFactory& factory,
        const AssemblyOptions& assemblyOptions)
{
    Fiber::SurfaceNormalDependentFunction<Functor> fiberFunction(functor);
    return gridFunctionFromFiberFunction(
                space, fiberFunction, factory, assemblyOptions);
}

} // namespace Bempp

#endif
