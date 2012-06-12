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


#ifndef bempp_interpolated_function_hpp
#define bempp_interpolated_function_hpp

#include "../common/common.hpp"

#include "../grid/vtk_writer.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../fiber/function.hpp"

#include "../common/armadillo_fwd.hpp"

namespace Bempp
{

class Grid;

template <typename ValueType> class InterpolatedFunction;

/** \brief Function defined by its values at a set of interpolation points
      and an interpolation method. */
template <typename ValueType>
class InterpolatedFunction : public Fiber::Function<ValueType>
{
public:
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

    enum InterpolationMethod {
        LINEAR
    };

    /** \brief Construct function given its values at vertices of a grid. */
    InterpolatedFunction(const Grid& grid,
                         const arma::Mat<ValueType>& vertexValues,
                         InterpolationMethod method = LINEAR);

    /** \brief Interpolation grid. */
    const Grid& grid() const;

    virtual size_t worldDimension() const;
    virtual size_t codomainDimension() const;
    virtual void addGeometricalDependencies(size_t& geomDeps) const;

    virtual void evaluate(const Fiber::GeometricalData<CoordinateType>& geomData,
                          arma::Mat<ValueType>& result) const;

//    virtual void evaluate(const arma::Mat<ValueType>& global,
//                          arma::Mat<ValueType>& values) const;

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
    void exportToVtk(const char* dataLabel,
                     const char* fileNamesBase, const char* filesPath = 0,
                     VtkWriter::OutputType type = VtkWriter::ASCII) const;

//    /** \brief Copy vertex values from a function defined on a subset of the
//      surface of the interpolation grid. */
//    void setSurfaceValues(const GridFunction<ValueType>& surfaceFunction);

//    /** \brief Copy vertex values from a function interpolated on a surface grid. */
//    void setSurfaceValues(const InterpolatedFunction<ValueType>& surfaceFunction);

    InterpolatedFunction<ValueType>& operator+=(
            const InterpolatedFunction<ValueType> &rhs);
    InterpolatedFunction<ValueType>& operator-=(
            const InterpolatedFunction<ValueType> &rhs);
    InterpolatedFunction<ValueType>& operator*=(ValueType rhs);
    InterpolatedFunction<ValueType>& operator/=(ValueType rhs);

    const InterpolatedFunction<ValueType> operator+(
            const InterpolatedFunction<ValueType> &other) const;
    const InterpolatedFunction<ValueType> operator-(
            const InterpolatedFunction<ValueType> &other) const;

    const InterpolatedFunction<ValueType> operator/(
            ValueType other) const;

private:
    void checkCompatibility(const InterpolatedFunction<ValueType>& other) const;

private:
    const Grid& m_grid;
    arma::Mat<ValueType> m_vertexValues;
    InterpolationMethod m_method;
};

template <typename ValueType>
const InterpolatedFunction<ValueType> operator*(
        ValueType lhs, const InterpolatedFunction<ValueType>& rhs)
{
    return InterpolatedFunction<ValueType>(rhs) *= lhs;
}

template <typename ValueType>
const InterpolatedFunction<ValueType> operator*(
        const InterpolatedFunction<ValueType>& lhs, ValueType rhs)
{
    return operator*(rhs, lhs);
}

} // namespace Bempp

#endif
