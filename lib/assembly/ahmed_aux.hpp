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

#include <memory>
#include <stdexcept>

// Ahmed's include files
#include <bemblcluster.h>
#include <matrix.h>
#include <matgen_sqntl.h>
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
        switch (dim)
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

template <typename ValueType, typename GeometryTypeRows, typename GeometryTypeCols>
class AhmedMatrix : public Matrix<ValueType>
{
public:
    typedef bemblcluster<GeometryTypeRows, GeometryTypeCols>
    AhmedBemblcluster;

    AhmedMatrix(unsigned int rowCount, unsigned int columnCount,
                std::auto_ptr<AhmedBemblcluster> blockCluster,
                mblock<ValueType>** blocks) :
        Matrix<double>(rowCount, columnCount),
        m_blockCluster(blockCluster), m_blocks(blocks)
    {
    }

    virtual ~AhmedMatrix() {
        freembls(m_blockCluster.get(), m_blocks);
    }

    void amux(ValueType d, ValueType* x, ValueType* y) const {
        multaHvec(d, m_blockCluster.get(), m_blocks, x, y);
    }

    void precond_apply(double* x) const {
        // TODO
    }

private:
    std::auto_ptr<AhmedBemblcluster> m_blockCluster;
    mblock<ValueType>** m_blocks;
};

} // namespace Bempp

#endif
