#ifndef bempp_grid_function_hpp
#define bempp_grid_function_hpp

#include "../grid/vtk_writer.hpp"

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

// OUT-OF-DATE COMMENT
/*
The GridFunction object has two envisaged applications.

First, its assembleWeakForm() method can be called to evaluate the projections
of basis functions of a given space on the grid function. This is used in the
solution of BIEs via the Galerkin method.

This application will essentially be handled by the
GridFunctionFromExactFunction subclass.

Second, it can be passed to the applyOffSurface() -- or, in future,
applyOnSurface() -- method of a linear operator. The operator will call the
GridFunction::evaluate() method. In this application, it is conceivable to use
both the GridFunctionFromExactFunction and GridFunctionFromCoefficients
subclasses (e.g. one for known Dirichlet data and the other for calculated
Neumann data). That's why it is important to have the same interface for both
subclasses.

Note that in future we will want to handle mixed problems (part-Dirichlet,
part-Neumann), and then subclasses implementing a combination of exact and
coefficient-derived functions will be needed. However, the interface will stay
the same.
*/

// TODO: We can remove "nonordinary" fiber::function.

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

    /** \brief Grid on which this function is defined. */
    const Grid& grid() const;

    /** \brief Space in which this function is expanded. */
    const Space<ValueType>& space() const;

    int codomainDimension() const;

    // possibly replace output type with DiscreteFunction/GridFunctionCoefficients/sth like this
    arma::Col<ValueType> coefficients() const;
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

} // namespace Bempp

#endif
