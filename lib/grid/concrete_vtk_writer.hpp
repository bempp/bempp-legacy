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

#ifndef bempp_concrete_vtk_writer_hpp
#define bempp_concrete_vtk_writer_hpp

#include "../common/common.hpp"

#include "vtk_writer.hpp"

#include "../common/armadillo_fwd.hpp"
#include <memory>
#include <string>

namespace Bempp
{

// Forward declarations
template<typename DuneGridView> class ConcreteGridView;

/** \brief Wrapper of a Dune VTK writer for a grid view of type \p DuneGridView. */
template <typename DuneGridView>
class ConcreteVtkWriter : public VtkWriter
{
    Dune::VTKWriter<DuneGridView> m_dune_vtk_writer;
    const DuneGridView* m_dune_gv;

    friend std::auto_ptr<VtkWriter> ConcreteGridView<DuneGridView>::vtkWriter(
            Dune::VTK::DataMode dm) const;

    /** \brief Construct a VtkWriter working on a specific \p DuneGridView.
     *
     *  \param dune_gv The grid view the grid functions live on.
     *    (E.g. a \p LevelGridView.)
     *  \param dm The data mode.
     *
     *  \internal This constructor can only be called by the factory method
     *  ConcreteGridView::vtkWriter().
     */
    explicit ConcreteVtkWriter(const DuneGridView &dune_gv,
                               Dune::VTK::DataMode dm=Dune::VTK::conforming) :
        m_dune_vtk_writer(dune_gv, dm), m_dune_gv(&dune_gv) {
    }

public:
    virtual void clear() {
        m_dune_vtk_writer.clear();
    }

    virtual std::string write(const std::string &name,
                              OutputType type = ASCII) {
        return m_dune_vtk_writer.write(name, duneVtkOutputType(type));
    }

    virtual std::string pwrite(const std::string& name,
                               const std::string& path,
                               const std::string& extendpath,
                               OutputType type = ASCII) {
        return m_dune_vtk_writer.pwrite(name, path, extendpath,
                                        duneVtkOutputType(type));
    }

private:
    virtual void addCellDataDoubleImpl(const arma::Mat<double>& data,
                                       const std::string &name) {
        addCellDataImpl(data, name);
    }

    virtual void addCellDataFloatImpl(const arma::Mat<float>& data,
                                      const std::string &name) {
        addCellDataImpl(data, name);
    }

    template <typename ValueType>
    void addCellDataImpl(const arma::Mat<ValueType>& data,
                              const std::string &name) {
        const size_t ncomp = data.n_rows;
        if (ncomp < 1)
            return; // empty matrix
        if ((int)data.n_cols != m_dune_gv->size(0 /* cell codim */))
            throw std::logic_error("VtkWriter::addCellData(): number of columns "
                                   "of 'data' different from the number of cells");
        m_dune_vtk_writer.addCellData(data, name, ncomp);
    }

    virtual void addVertexDataDoubleImpl(const arma::Mat<double>& data,
                                       const std::string &name) {
        addVertexDataImpl(data, name);
    }

    virtual void addVertexDataFloatImpl(const arma::Mat<float>& data,
                                      const std::string &name) {
        addVertexDataImpl(data, name);
    }

    template <typename ValueType>
    void addVertexDataImpl(const arma::Mat<ValueType>& data,
                           const std::string &name) {
        const size_t ncomp = data.n_rows;
        if (ncomp < 1)
            return; // empty matrix
        if ((int)data.n_cols != m_dune_gv->size(DuneGridView::dimension /* vertex codim */))
            throw std::logic_error("VtkWriter::addVertexData(): number of columns "
                                   "of 'data' different from the number of vertices");
        m_dune_vtk_writer.addVertexData(data, name, ncomp);
    }

    Dune::VTK::OutputType duneVtkOutputType(OutputType type) const
    {
        switch (type)
        {
        case ASCII:
            return Dune::VTK::ascii;
        case BASE_64:
            return Dune::VTK::base64;
        case APPENDED_RAW:
            return Dune::VTK::appendedraw;
        case APPENDED_BASE_64:
            return Dune::VTK::appendedbase64;
        default:
            return static_cast<Dune::VTK::OutputType>(type);
        }
    }
};

} // namespace Bempp

#endif
