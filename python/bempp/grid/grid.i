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


		 %}


}

} // namespace Bempp

%include "grid/grid.hpp"
