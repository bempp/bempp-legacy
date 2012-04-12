#ifndef bempp_interpolated_function_hpp
#define bempp_interpolated_function_hpp

#include "../grid/vtk_writer.hpp"
#include "../fiber/function.hpp"

#include <armadillo>

namespace Bempp
{

class Grid;
template <typename ValueType> class GridFunction;

template <typename ValueType> class InterpolatedFunction;

template <typename ValueType>
const InterpolatedFunction<ValueType> operator*(
        ValueType lhs, const InterpolatedFunction<ValueType>& rhs);

/** \brief Function defined by its values at a set of interpolation points
      and an interpolation method. */
template <typename ValueType>
class InterpolatedFunction : public Fiber::Function<ValueType>
{
public:
    enum InterpolationMethod {
        LINEAR
    };

    /** \brief Construct function given its values at vertices of a grid. */
    InterpolatedFunction(const Grid& grid,
                         const arma::Mat<ValueType>& vertexValues,
                         InterpolationMethod method = LINEAR);

    /** \brief Interpolation grid. */
    const Grid& grid() const;

    virtual int worldDimension() const;
    virtual int codomainDimension() const;
    virtual void addGeometricalDependencies(int& geomDeps) const;

    virtual void evaluate(const Fiber::GeometricalData<ValueType>& geomData,
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

    /** \brief Copy vertex values from a function defined on a subset of the
      surface of the interpolation grid. */
    void setSurfaceValues(const GridFunction<ValueType>& surfaceFunction);

    /** \brief Copy vertex values from a function interpolated on a surface grid. */
    void setSurfaceValues(const InterpolatedFunction<ValueType>& surfaceFunction);

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

    const InterpolatedFunction<ValueType> operator*(
            ValueType other) const;
    const InterpolatedFunction<ValueType> operator/(
            ValueType other) const;

//    template <typename ValueType>
    friend const InterpolatedFunction<ValueType> operator*(
            ValueType lhs, const InterpolatedFunction<ValueType>& rhs);

private:
    void checkCompatibility(const InterpolatedFunction<ValueType>& other) const;

private:
    const Grid& m_grid;
    arma::Mat<ValueType> m_vertexValues;
    InterpolationMethod m_method;
};

} // namespace Bempp

#endif
