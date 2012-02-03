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

#ifndef bempp_geometry_adapter_hpp
#define bempp_geometry_adapter_hpp

#include "geometry.hpp"
#include "../fiber/geometry.hpp"

namespace Bempp
{

/** \brief Wrapper of Bempp::Geometry conforming to the Fiber::Geometry interface. */
class GeometryAdapter : public Fiber::Geometry<ctype, GeometryAdapter>
{
public:
    explicit GeometryAdapter(const Geometry& geometry) :
        m_wrapped(geometry) {
    }

    int dimension() const  {
        return asImp().dim();
    }

    void getInformation(
            int pointCount,
            const ctype* localCoords,
            ctype* globalCoords,
            ctype* normals,
            ctype* integrationElements,
            ctype* jacobiansTransposed,
            ctype* jacobianInversesTransposed) const {
        const int worldDim = worldDimension();
        const int dim = dimension();

        // The const_cast is there so that the constant matrices passed to
        // Bempp::Geometry member functions can be declared in the usual way,
        // as "const arma::Mat<ctype>&", instead of the clumsier variant
        // "const arma::Mat<const ctype>&".
        const arma::Mat<ctype> arma_localCoords(
                    const_cast<ctype*>(localCoords), dim, pointCount,
                    false /* copy_mem */);
        if (globalCoords)
        {
            arma::Mat<ctype> arma_globalCoords(
                        globalCoords, worldDim, pointCount, false);
            asImp().local2global(arma_localCoords, arma_globalCoords);
        }

        if (jacobiansTransposed)
        {
            arma::Cube<ctype> arma_jacobiansTransposed(
                        jacobiansTransposed, dim, worldDim, pointCount, false);
            asImp().jacobianTransposed(arma_localCoords, arma_jacobiansTransposed);
            if (normals)
                ; // TODO: calculate normals
        }
        else if (normals)
        {
            arma::Cube<ctype> arma_jacobiansTransposed(
                        dim, worldDim, pointCount);
            asImp().jacobianTransposed(arma_localCoords, arma_jacobiansTransposed);
            // TODO: calculate normals
        }

        if (jacobianInversesTransposed)
        {
            arma::Cube<ctype> arma_jacobianInversesTransposed(
                        jacobianInversesTransposed, worldDim, dim, pointCount, false);
            asImp().jacobianInverseTransposed(arma_localCoords,
                                              arma_jacobianInversesTransposed);
        }

        if (integrationElements)
        {
            arma::Row<ctype> arma_integrationElements(
                        integrationElements, pointCount, false);
            asImp().integrationElement(arma_localCoords, arma_integrationElements);
        }
    }

    void normal(const arma::Mat<ctype>& local,
                arma::Mat<ctype>& normals) const;
private:
    const Geometry& asImp() const {
        return m_wrapped;
    }

    const Geometry& m_wrapped;
};

} // namespace Bempp

#endif
