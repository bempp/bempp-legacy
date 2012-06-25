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
        self.vertices={}
        for v in vertex_view:
            id=id_set.entityId(v)
            self.vertices[id]=v.geometry().corners()[:,0]

        element_view=view.entities(0)
        self.elements={}
        for e in element_view:
            points=e.subEntities(2)
            elem_id=id_set.entityId(e)
            self.elements[elem_id]=[]
            for p in points:
                self.elements[elem_id].append(id_set.entityId(p))

    def getVertices(self):
        """Return the vertices of the grid"""
        if not hasattr(self,'vertices'): self.__enumerateGridData()
        return self.vertices

    def getElements(self):
        """Return the elements of the grid"""
        if not hasattr(self,'elements'): self.__enumerateGridData()
        return self.elements

    def plot(self):
        """Plot the grid using VTK"""
        try:
            import vtk
        except ImportError:
            print "The Python VTK bindings needs to be installed for this method!"
            return

        if not hasattr(self,'vertices'): self.__enumerateGridData()

        def create_cell(elem):
            polygon=vtk.vtkPolygon()
            ids=polygon.GetPointIds()
            ids.SetNumberOfIds(len(elem))
            for i,p in enumerate(elem):
                ids.SetId(i,p)
            return polygon

        points=vtk.vtkPoints()
        points.SetNumberOfPoints(len(self.vertices.keys()))
        for i in self.vertices:
            points.InsertPoint(i,self.vertices[i][0],self.vertices[i][1],self.vertices[i][2])
        polyGrid=vtk.vtkUnstructuredGrid()
        polyGrid.SetPoints(points)
        cells=[create_cell(self.elements[i]) for i in self.elements]
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
