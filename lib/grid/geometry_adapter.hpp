#ifndef bempp_geometry_adapter_hpp
#define bempp_geometry_adapter_hpp

#include "index_set.hpp"
#include "entity.hpp"
#include "geometry.hpp"

#include <vector>

/** Wrapper around Geometry conforming to the interface required by Fiber. */

namespace Bempp
{

class GeometryAdapter
{
public:
    typedef IndexSet::IndexType IndexType;

    GeometryAdapter(const Entity<0>& element, const IndexSet& indexSet) :
        m_geometry(element.geometry()) {
        const int codimVertex = m_geometry.dim();
        const int vertexCount = m_geometry.cornerCount();
        m_vertexIndices.resize(vertexCount);
        // To optimise
        for (int i = 0; i < vertexCount; ++i)
            m_vertexIndices[i] = indexSet.subEntityIndex(element, codimVertex, i);
    }

    int dimension() const {
        return m_geometry.dim();
    }

    IndexType vertexIndex(int i) const {
        return m_vertexIndices[i];
    }

    int vertexCount() const {
        return m_vertexIndices.size();
    }

    void getData(int what, const arma::Mat<ctype>& local,
                 Fiber::GeometricalData<ctype>& data) const {
        m_geometry.getData(what, local, data);
    }

private:
    const Geometry& m_geometry;
    std::vector<IndexType> m_vertexIndices;
};

} // namespace Bempp

#endif
