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

#include "../grid/vtk_writer.hpp"
#include "../fiber/scalar_traits.hpp"
#include "vector.hpp"

#include <armadillo>
#include <memory>

namespace Fiber
{

template <typename BasisValueType, typename ResultType, typename GeometryFactory>
class LocalAssemblerFactory;
template <typename ValueType> class Basis;
template <typename ValueType> class Function;
template <typename ValueType> class LocalAssemblerForGridFunctions;

} // namespace Fiber

namespace Bempp
{

class AssemblyOptions;
class GeometryFactory;
class Grid;
template <int codim> class Entity;
template <typename ValueType> class Space;

/** \brief Function defined on a grid. */
template <typename BasisFunctionType, typename ResultType>
class GridFunction
{
public:
    typedef Fiber::LocalAssemblerFactory<BasisFunctionType, ResultType, GeometryFactory>
    LocalAssemblerFactory;
    typedef Fiber::LocalAssemblerForGridFunctions<ResultType> LocalAssembler;
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

    /** \brief Construct by evaluating the expansion coefficients of a global
      function in the provided function space. */
    GridFunction(const Space<BasisFunctionType>& space,
                 const Fiber::Function<ResultType>& function,
                 const LocalAssemblerFactory& factory,
                 const AssemblyOptions& assemblyOptions);

    /** \brief Construct from known expansion coefficients in the provided function space. */
    GridFunction(const Space<BasisFunctionType>& space,
                 const arma::Col<ResultType>& coefficients);

    /** \brief Construct from known expansion coefficients in the provided function space,
        represented as Vector<ValueType>. */
    GridFunction(const Space<BasisFunctionType>& space,
                 const Vector<ResultType>& coefficients);

    /** \brief Grid on which this function is defined. */
    const Grid& grid() const;

    /** \brief Space in which this function is expanded. */
    const Space<BasisFunctionType>& space() const;

    int codomainDimension() const;

    Vector<ResultType> coefficients() const;
    void setCoefficients(const Vector<ResultType>& coeffs);

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
    /** \brief Calculate projections of the function on test functions from
      the given space. */
    arma::Col<ResultType> calculateProjections(
            const Fiber::Function<ResultType>& globalFunction,
            const Space<BasisFunctionType>& space,
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

    arma::Col<ResultType> reallyCalculateProjections(
            const Space<BasisFunctionType>& space,
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    /** \brief Evaluate function at either vertices or barycentres. */
    void evaluateAtSpecialPoints(
            VtkWriter::DataType dataType, arma::Mat<ResultType>& result) const;

private:
    const Space<BasisFunctionType>& m_space;
    arma::Col<ResultType> m_coefficients;
};

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

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const ScalarType& scalar, const GridFunction<BasisFunctionType, ResultType>& g2);

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator/(
        const GridFunction<BasisFunctionType, ResultType>& g1, const ScalarType& scalar);

} // namespace Bempp

#endif
