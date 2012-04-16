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
#include "vector.hpp"

#include <armadillo>
#include <memory>

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
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
template <typename ValueType>
class GridFunction
{
public:
    typedef Fiber::LocalAssemblerFactory<ValueType, GeometryFactory>
    LocalAssemblerFactory;
    typedef Fiber::LocalAssemblerForGridFunctions<ValueType> LocalAssembler;

    /** \brief Construct by evaluating the expansion coefficients of a global
      function in the provided function space. */
    GridFunction(const Space<ValueType>& space,
                 const Fiber::Function<ValueType>& function,
                 const LocalAssemblerFactory& factory,
                 const AssemblyOptions& assemblyOptions);

    /** \brief Construct from known expansion coefficients in the provided function space. */
    GridFunction(const Space<ValueType>& space,
                 const arma::Col<ValueType>& coefficients);

    /** \brief Construct from known expansion coefficients in the provided function space,
        represented as Vector<ValueType>. */
    GridFunction(const Space<ValueType>& space,
                 const Vector<ValueType>& coefficients);

    /** \brief Grid on which this function is defined. */
    const Grid& grid() const;

    /** \brief Space in which this function is expanded. */
    const Space<ValueType>& space() const;

    int codomainDimension() const;

    // possibly replace output type with DiscreteFunction/GridFunctionCoefficients/sth like this
    Vector<ValueType> coefficients() const;
    const Fiber::Basis<ValueType>& basis(const Entity<0>& element) const;
    void getLocalCoefficients(const Entity<0>& element,
                              std::vector<ValueType>& coeffs) const;

    /** Export the function to a VTK file.

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
    arma::Col<ValueType> calculateProjections(
            const Fiber::Function<ValueType>& globalFunction,
            const Space<ValueType>& space,
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

    arma::Col<ValueType> reallyCalculateProjections(
            const Space<ValueType>& space,
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    /** \brief Evaluate function at either vertices or barycentres. */
    void evaluateAtSpecialPoints(
            VtkWriter::DataType dataType, arma::Mat<ValueType>& result) const;

private:
    const Space<ValueType>& m_space;
    arma::Col<ValueType> m_coefficients;
};

template<typename ValueType>
GridFunction<ValueType> operator+(const GridFunction<ValueType>& g1, const GridFunction<ValueType>& g2);

template<typename ValueType>
GridFunction<ValueType> operator-(const GridFunction<ValueType>& g1, const GridFunction<ValueType>& g2);

template<typename ValueType>
GridFunction<ValueType> operator*(const GridFunction<ValueType>& g1, const ValueType& scalar);

template<typename ValueType>
GridFunction<ValueType> operator*(const ValueType& scalar, const GridFunction<ValueType>& g2);

template<typename ValueType>
GridFunction<ValueType> operator/(const GridFunction<ValueType>& g1, const ValueType& scalar);

} // namespace Bempp

#endif
