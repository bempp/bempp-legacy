%{
#include "grid/grid.hpp"
#include "grid/grid_view.hpp"
#include "grid/entity_iterator.hpp"
%}

%include "grid_docstrings.i"

namespace Bempp 
{
    
%extend Grid 
{
    %pythonappend leafView %{
        val._parentGrid = self
    %}

    %pythonappend levelView %{
        val._parentGrid = self
    %}

    %pythonappend globalIdSet %{
        val._parentGrid = self
    %}

    // this function is only for internal use
    %ignore elementGeometryFactory;

   %pythoncode %{
    def __enumerateGridData(self):
        view=self.leafView()
        id_set=self.globalIdSet()
        vertex_view=view.entities(2)
        self.__vertices={}
        for v in vertex_view:
            id=id_set.entityId(v)
            self.__vertices[id]=v.geometry().corners()[:,0]

        element_view=view.entities(0)
        self.__elements={}
        for e in element_view:
            points=e.subEntities(2)
            elem_id=id_set.entityId(e)
            self.__elements[elem_id]=[]
            for p in points:
                self.__elements[elem_id].append(id_set.entityId(p))

    def getVertices(self):
        """Return the vertices of the grid"""
        if not hasattr(self,'_Grid__vertices'): self.__enumerateGridData()
        return self.__vertices

    def getElements(self):
        """Return the elements of the grid"""
        if not hasattr(self,'_Grid__elements'): self.__enumerateGridData()
        return self.__elements

    def plot(self):
        """Plot the grid using VTK"""
        try:
            import vtk
        except ImportError:
            print "The Python VTK bindings needs to be installed for this method!"
            return

        if not hasattr(self,'__Grid_vertices'): self.__enumerateGridData()

        def create_cell(elem):
            polygon=vtk.vtkPolygon()
            ids=polygon.GetPointIds()
            ids.SetNumberOfIds(len(elem))
            for i,p in enumerate(elem):
                ids.SetId(i,p)
            return polygon

        points=vtk.vtkPoints()
        points.SetNumberOfPoints(len(self.__vertices.keys()))
        for i in self.__vertices:
            points.InsertPoint(i,self.__vertices[i][0],self.__vertices[i][1],self.__vertices[i][2])
        polyGrid=vtk.vtkUnstructuredGrid()
        polyGrid.SetPoints(points)
        cells=[create_cell(self.__elements[i]) for i in self.__elements]
        for cell in cells:
            polyGrid.InsertNextCell(cell.GetCellType(),cell.GetPointIds())
        mapper=vtk.vtkDataSetMapper()
        mapper.SetInput(polyGrid)
        actor=vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetRepresentationToWireframe()
        renderer=vtk.vtkRenderer()
        renderer.AddActor(actor)
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

} // namespace Bempp

%include "grid/grid.hpp"
