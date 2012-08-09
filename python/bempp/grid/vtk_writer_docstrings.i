// Docstrings ------------------------------------------------------------------

%define VtkWriter_docstring
"Exporter of data in the vtk format.

Exports data (living on cells or vertices of a grid) to a file
suitable for easy visualization with the Visualization Toolkit (VTK;
see http://public.kitware.com/VTK)"
%enddef 

%define VtkWriter_addCellData_docstring
"Add to the visualization output a grid function representing data
associated with the cells of the grid.

*Arguments*
    data (ndarray)
        2D array whose (m, n)th entry contains the value of the mth component
        of the grid function in the nth cell.
    name (string)
        Name to identify the grid function."
%enddef

%define VtkWriter_addVertexData_docstring
"Add to the visualization output a grid function representing data
associated with the vertices of the grid.

*Arguments*
    data (ndarray)
        2D array whose (m, n)th entry contains the value of the mth component
        of the grid function at the nth vertex.
    name (string)
        Name to identify the grid function."
%enddef

%define VtkWriter_clear_docstring
"Clear the list of registered functions."
%enddef

%define VtkWriter_write_autodoc_docstring
"write(self, name, type = 'ascii') -> string"
%enddef

%define VtkWriter_write_docstring
"Write output (interface might change later).

This method can be used in parallel as well as in serial programs.
For serial runs (commSize=1) it chooses other names without the
`s####:p####:' prefix for the .vtu/.vtp files and omits writing of the
.pvtu/pvtp file however.  For parallel runs (commSize > 1) it is the
same as a call to pwrite() with path='' and extendpath=''.

*Arguments:*
    name (string)
        Basic name to write (may not contain a path).
    type (string, optional)
        Format of the output, one of `ascii' (default), `base64',
        `appendedraw' or `appendedbase64'.
 
Returns the name of the created file."
%enddef

%define VtkWriter_pwrite_autodoc_docstring
"pwrite(self, name, path, extendpath, type = 'ascii') -> string"
%enddef

%define VtkWriter_pwrite_docstring
"Write output (interface might change later).

`pwrite' means `path write' (i.e. write somewhere else than the current
directory). The `p' does not mean this method has a monopoly on
parallel writing, the regular write() method can do that just fine.

*Arguments:*
    name (string)        
        Base name of the output files. This should not contain any
        directory part or filename extensions. It will be used
        both for the piece file of each process and the parallel
        collection file.
    path (string)        
        Directory where to put the parallel collection (.pvtu/.pvtp)
        file. If it is relative, it is taken relative to the current
        directory.
    extendpath (string)                
        Directory where to put the piece file (.vtu/.vtp) of this
        process. If it is relative, it is taken relative to the
        directory denoted by path.
    type (string, optional)
        Format of the output, one of `ascii' (default), `base64',
        `appendedraw' or `appendedbase64'.
 
Returns the name of the created file.

*Note:* Currently, extendpath may not be absolute unless path is."
%enddef

// Declarations ----------------------------------------------------------------

namespace Bempp
{

DECLARE_CLASS_DOCSTRING (VtkWriter);
DECLARE_METHOD_DOCSTRING(VtkWriter, addCellData, 1);
DECLARE_METHOD_DOCSTRING(VtkWriter, addVertexData, 1);
DECLARE_METHOD_DOCSTRING(VtkWriter, clear, 1);
DECLARE_METHOD_DOCSTRING(VtkWriter, write, 0);
DECLARE_METHOD_DOCSTRING(VtkWriter, pwrite, 0);

} // namespace Bempp
