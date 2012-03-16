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

#ifndef bempp_ahmed_aux_hpp
#define bempp_ahmed_aux_hpp

#include "../common/types.hpp"

#include <boost/scoped_array.hpp>
#include <boost/shared_array.hpp>
#include <iostream>
#include <memory>
#include <stdexcept>

// Ahmed's include files
#include <apprx.h>
#include <bemblcluster.h>
#include <bllist.h>
// #include <matgen_sqntl.h>
// #include <matgen_omp.h>
#include <matrix.h>
#undef SIGN
#undef SQR
#undef MIN
#undef MAX

namespace Bempp {

template <typename ValueType>
struct AhmedDofWrapper : public Point3D<ValueType>
{
    ValueType getcenter(int dim) const
    {
        switch  (dim)
        {
        case 0: return this->x;
        case 1: return this->y;
        case 2: return this->z;
        default:
            throw std::invalid_argument("AhmedDofWrapper::getcenter(): "
                                        "invalid dimension index");
        }
    }
};

class AhmedMblockArrayDeleter
{
public:
    AhmedMblockArrayDeleter(size_t arraySize) :
        m_arraySize(arraySize) {
    }

    template <typename ValueType>
    void operator() (mblock<ValueType>** blocks) const {
        freembls(m_arraySize, blocks);
    }

private:
    size_t m_arraySize;
};

template <typename ValueType>
boost::shared_array<mblock<ValueType>*> allocateAhmedMblockArray(
        const blcluster* cluster)
{
    mblock<ValueType>** blocks = 0;
    const size_t blockCount = cluster->nleaves();
    allocmbls(blockCount, blocks);
    return boost::shared_array<mblock<ValueType>*>(
                blocks, AhmedMblockArrayDeleter(blockCount));
}

class AhmedLeafClusterArray
{
public:
    // The parameter should actually be const, but Ahmed's gen_BlSequence lacks
    // const-correctness
    explicit AhmedLeafClusterArray(blcluster* clusterTree) :
        m_size(0) {
        blcluster** leafClusters = 0;
        try {
            gen_BlSequence(clusterTree, leafClusters);
        }
        catch (...) {
            delete[] leafClusters;
            throw; // rethrow the exception
        }
        m_leafClusters.reset(leafClusters);
        m_size = clusterTree->nleaves();
    }

    size_t size() const {
        return m_size;
    }

    blcluster* operator[] (size_t n) {
        return m_leafClusters[n];
    }

    const blcluster* operator[] (size_t n) const {
        return m_leafClusters[n];
    }

private:
    boost::scoped_array<blcluster*> m_leafClusters;
    size_t m_size;
};

template <typename ValueType, typename GeometryTypeRows, typename GeometryTypeCols>
class AhmedMatrix : public Matrix<ValueType>
{
public:
    typedef bemblcluster<GeometryTypeRows, GeometryTypeCols>
    AhmedBemblcluster;

    AhmedMatrix(unsigned int rowCount, unsigned int columnCount,
                std::auto_ptr<AhmedBemblcluster> blockCluster,
                boost::shared_array<mblock<ValueType>*> blocks) :
        Matrix<double>(rowCount, columnCount),
        m_blockCluster(blockCluster), m_blocks(blocks)
    {
    }

    virtual void amux(ValueType d, ValueType* x, ValueType* y) const {
        AhmedBemblcluster* ptr = m_blockCluster.get();
        multaHvec(d, m_blockCluster.get(), m_blocks.get(), x, y);
    }

    virtual void precond_apply(double* x) const {
        // TODO
    }

private:
    std::auto_ptr<AhmedBemblcluster> m_blockCluster;
    boost::shared_array<mblock<ValueType>*> m_blocks;
};

} // namespace Bempp

#endif
