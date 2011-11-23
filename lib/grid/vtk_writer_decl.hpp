// Copyright (C) 2011 by the BEM++ Authors
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

#ifndef bempp_vtk_writer_decl_hpp
#define bempp_vtk_writer_decl_hpp

#include <armadillo>
#include <memory>
#include <string>

namespace Bempp
{

/** \brief Writer for the output of data in the vtk format.
  *
  * Writes data (living on cells or vertices of a grid)
  * to a file suitable for easy visualization with
  * <a href="http://public.kitware.com/VTK/">The Visualization Toolkit (VTK)</a>.
  */
class VtkWriter
{
public:
    /** Destructor */
    virtual ~VtkWriter() {}

    /** \brief Add a grid function (represented by container) that lives on the cells of
     *  the grid to the visualization output.
     *
     *  \param data A matrix whose (m,n)th entry contains the value of the mth component of the grid function in the nth cell.
     *  \param name A name to identify the grid function.
     */
    virtual void addCellData(const arma::Mat<double>& data, const std::string &name) = 0;

    /** \brief Add a grid function (represented by container) that lives on the vertices of the
     *  grid to the visualization output.
     *
     *  \param data A matrix whose (m,n)th entry contains the value of the mth component of the grid function in the nth vertex.
     *  \param name A name to identify the grid function.
     */
    virtual void addVertexData(const arma::Mat<double>& data, const std::string &name) = 0;

    /** Clear the list of registered functions. */
    virtual void clear() = 0;

    /** \brief write output (interface might change later)
     *
     *  This method can be used in parallel as well as in serial programs.
     *  For serial runs (commSize=1) it chooses other names without the
     *  "s####:p####:" prefix for the .vtu/.vtp files and omits writing of the
     *  .pvtu/pvtp file however.  For parallel runs (commSize > 1) it is the
     *  same as a call to pwrite() with path="" and extendpath="".
     *
     *  \param[in]  name  Basic name to write (may not contain a path).
     *  \param[in]  type  Type of output (e.g,, ASCII) (optional).
     *
     *  \returns Name of the created file.
     */
    virtual std::string write (const std::string &name,
                               Dune::VTK::OutputType type = Dune::VTK::ascii) = 0;

    /** \brief write output (interface might change later)
     *
     * "pwrite" means "path write" (i.e. write somewhere else than the current
     * directory).  The "p" does not mean this method has a monopoly on
     * parallel writing, the regular write(const std::string &,
     * VTK::OutputType) method can do that just fine.
     *
     * \param name       Base name of the output files.  This should not
     *                   contain any directory part and not filename
     *                   extensions.  It will be used both for each processes
     *                   piece as well as the parallel collection file.
     * \param path       Directory where to put the parallel collection
     *                   (.pvtu/.pvtp) file.  If it is relative, it is taken
     *                   realtive to the current directory.
     * \param extendpath Directory where to put the piece file (.vtu/.vtp) of
     *                   this process.  If it is relative, it is taken
     *                   relative to the directory denoted by path.
     * \param type       How to encode the data in the file.
     *
     * \returns Name of the created file.
     *
     * \note Currently, extendpath may not be absolute unless path is
     *       absolute, because that would require the value of the current
     *       directory.
     *
     * \throw NotImplemented Extendpath is absolute but path is relative.
     * \throw IOError        Failed to open a file.
     */
    virtual std::string pwrite(const std::string& name, const std::string& path, const std::string& extendpath,
                               Dune::VTK::OutputType type = Dune::VTK::ascii) = 0;
};

template <typename DuneGridView>
class ConcreteVtkWriter : public VtkWriter
{
    Dune::VTKWriter<DuneGridView> m_dune_vtk_writer;
    const DuneGridView* m_dune_gv;

    friend std::auto_ptr<VtkWriter> ConcreteGridView<DuneGridView>::vtkWriter(Dune::VTK::DataMode dm) const;

    /** \brief Construct a VtkWriter working on a specific DuneGridView.
     *
     *  \param gv The grid view the grid functions live on. (E.g. a LevelGridView.)
     *  \param dm The data mode.
     *
     *  \internal This constructor can only be called by the factory method ConcreteGridView::vtkWriter().
     */
    explicit ConcreteVtkWriter(const DuneGridView &dune_gv, Dune::VTK::DataMode dm=Dune::VTK::conforming) :
        m_dune_vtk_writer(dune_gv, dm), m_dune_gv(&dune_gv) {
    }

public:
    virtual void addCellData(const arma::Mat<double>& data, const std::string &name) {
        const int ncomp = data.n_rows;
        if (ncomp < 1)
            return; // empty matrix
        if (data.n_cols != m_dune_gv->size(0 /* cell codim */))
            throw std::logic_error("VtkWriter::addCellData(): number of columns of 'data' different from the number of cells");
        m_dune_vtk_writer.addCellData(data, name, ncomp);
    }

    virtual void addVertexData(const arma::Mat<double>& data, const std::string &name) {
        const int ncomp = data.n_rows;
        if (ncomp < 1)
            return; // empty matrix
        if (data.n_cols != m_dune_gv->size(DuneGridView::dimension /* vertex codim */))
            throw std::logic_error("VtkWriter::addVertexData(): number of columns of 'data' different from the number of vertices");
        m_dune_vtk_writer.addVertexData(data, name, ncomp);
    }

    virtual void clear() {
        m_dune_vtk_writer.clear();
    }

    virtual std::string write(const std::string &name,
                              Dune::VTK::OutputType type = Dune::VTK::ascii) {
        return m_dune_vtk_writer.write(name, type);
    }

    virtual std::string pwrite(const std::string& name, const std::string& path, const std::string& extendpath,
                               Dune::VTK::OutputType type = Dune::VTK::ascii) {
        return m_dune_vtk_writer.pwrite(name, path, extendpath, type);
    }
};

} // namespace Bempp

#endif
