%{
#include "assembly/grid_function.hpp"
%}

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);

%apply arma::Mat<float>& ARGOUT_MAT {
arma::Mat<float>& result_
};

%apply arma::Mat<double>& ARGOUT_MAT {
arma::Mat<double>& result_
};

%apply arma::Mat<std::complex<float> >& ARGOUT_MAT {
arma::Mat<std::complex<float> >& result_
};

%apply arma::Mat<std::complex<double> >& ARGOUT_MAT {
arma::Mat<std::complex<double> >& result_
};

%extend GridFunction
{
    %ignore setCoefficients;
    %ignore setProjections;
    %ignore codomainDimension;
    %ignore space;

    %ignore basis;
    %ignore getLocalCoefficients;

    %apply arma::Col<float>& ARGOUT_COL {
        arma::Col<float>& col_out
    };
    %apply arma::Col<double>& ARGOUT_COL {
        arma::Col<double>& col_out
    };
    %apply arma::Col<std::complex<float> >& ARGOUT_COL {
        arma::Col<std::complex<float> >& col_out
    };
    %apply arma::Col<std::complex<double> >& ARGOUT_COL {
        arma::Col<std::complex<double> >& col_out
    };

    void coefficients(arma::Col<ResultType>& col_out)
    {
        col_out = $self->coefficients();
    }

    void projections(arma::Col<ResultType>& col_out)
    {
        col_out = $self->projections();
    }

    %ignore coefficients;
    %ignore projections;

    %pythoncode %{

def plot(self,mode=VtkWriter.VERTEX_DATA,realImag='real'):
    """Plot a grid function in a VTK Window"""

    import numpy

    try:
        import vtk
    except ImportError:
        print "The Python VTK bindings needs to be installed for this method!"
        return

    if not mode in [VtkWriter.VERTEX_DATA, VtkWriter.CELL_DATA]:
        raise ValueError('Unknown mode specified. Valid modes are bempp.VtkWriter.VERTEX_DATA and bempp.VtkWriter.CELL_DATA!')

    if not realImag in ['real', 'imag']:
        raise ValueError("Unknown value for 'realImag'. It needs to be either 'real' or 'imag'!")


    polyGrid = self.grid().getVtkGrid()
    if mode==VtkWriter.VERTEX_DATA:
        values = self.evaluateAtSpecialPoints(VtkWriter.VERTEX_DATA).flatten()
        vtk_data = polyGrid.GetPointData()
    elif mode==VtkWriter.CELL_DATA:
        values = self.evaluateAtSpecialPoints(VtkWriter.CELL_DATA).flatten()
        vtk_data = polyGrid.GetCellData()

    if realImag=='real':
        values = numpy.real(values)
    else:
        values = numpy.imag(values)

    count = len(values)
    vtk_array = vtk.vtkDoubleArray()
    vtk_array.SetNumberOfValues(count)
    for i,p in enumerate(values): vtk_array.SetValue(i,p)
    vtk_data.SetScalars(vtk_array)

    mapper = vtk.vtkDataSetMapper()
    mapper.SetInput(polyGrid)
    mapper.SetScalarRange(vtk_array.GetRange())

    data_actor = vtk.vtkActor()
    data_actor.SetMapper(mapper)

    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(mapper.GetLookupTable())

    renderer=vtk.vtkRenderer()
    renderer.AddActor(data_actor)
    renderer.AddActor(scalar_bar)
    renderer.SetBackground(.1,.2,.4)

    window=vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    window.SetSize(800,600)

    irenderer=vtk.vtkRenderWindowInteractor()
    irenderer.SetRenderWindow(window)
    irenderer.Initialize()

    window.Render()
    irenderer.Start()

%}

}

%ignore gridFunctionFromFiberFunction;

} // namespace Bempp

%include "assembly/grid_function.hpp";

namespace Bempp
{

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction)
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    gridFunctionFromSurfaceNormalIndependentFunctor);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    gridFunctionFromSurfaceNormalDependentFunctor);

%clear arma::Mat<float>& result_;
%clear arma::Mat<double>& result_;
%clear arma::Mat<std::complex<float> >& result_;
%clear arma::Mat<std::complex<double> >& result_;

%clear arma::Col<float>& col_out;
%clear arma::Col<double>& col_out;
%clear arma::Col<std::complex<float> >& col_out;
%clear arma::Col<std::complex<double> >& col_out;

} // namespace Bempp

%pythoncode %{

def gridFunctionFromSurfaceNormalDependentFunctor(space, functor, factory,
        assemblyOptions):
    basisFunctionType = checkType(space.basisFunctionType())
    resultType = checkType(functor.valueType())
    result = constructObjectTemplatedOnBasisAndResult(
        "gridFunctionFromSurfaceNormalDependentFunctor",
        basisFunctionType, resultType,
        space, functor, factory, assemblyOptions)
    result._space = space
    return result

def gridFunctionFromSurfaceNormalIndependentFunctor(space, functor, factory,
        assemblyOptions, basisFunctionType='float64'):
    basisFunctionType = checkType(space.basisFunctionType())
    resultType = checkType(functor.valueType())
    result = constructObjectTemplatedOnBasisAndResult(
        "gridFunctionFromSurfaceNormalIndependentFunctor",
        basisFunctionType, resultType,
        space, functor, factory, assemblyOptions)
    result._space = space
    return result

%}
