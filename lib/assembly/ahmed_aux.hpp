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

// to ensure there are no inconsistent forward declarations
#include "ahmed_aux_fwd.hpp"

#include "ahmed_leaf_cluster_array.hpp"
#include "../common/types.hpp"

#include <boost/scoped_array.hpp>
#include <boost/shared_array.hpp>
#include <complex>
#include <iostream>
#include <memory>
#include <stdexcept>

// Ahmed's include files
#include <cmplx.h>

inline comp<float> operator*=(comp<float>& a, double b)
{
    return operator*=(a, static_cast<float>(b));
}

inline comp<float> operator/(comp<float>& a, double b)
{
    return operator/(a, static_cast<float>(b));
}

inline comp<float> operator/(double a, comp<float>& b)
{
    return operator/(static_cast<float>(a), b);
}

#include <apprx.h>
#include <bemblcluster.h>
#include <bllist.h>
// #include <matgen_sqntl.h>
// #define _OPENMP
#include <matgen_omp.h>
// #undef _OPENMP

bool multaHvec_omp(double d, blcluster* bl, mblock<double>** A, double* x,
           double* y);

#include <matrix.h>
#undef SIGN
#undef SQR
#undef MIN
#undef MAX

namespace Bempp
{

// Casts.

inline float ahmedCast(float x) {
    return x;
}

inline double ahmedCast(double x) {
    return x;
}

inline scomp ahmedCast(std::complex<float> x) {
    return scomp(x.real(), x.imag());
}

inline dcomp ahmedCast(std::complex<double> x) {
    return dcomp(x.real(), x.imag());
}

/** \brief An Ahmed-compatible degree-of-freedom type. */
template <typename CoordinateType>
struct AhmedDofWrapper : public Point3D<CoordinateType>
{
    CoordinateType getcenter(int dim) const
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

template <typename T>
class ExtendedBemCluster : public bemcluster<T>
{
public:
    ExtendedBemCluster(
            T* dofs, unsigned int* op_perm,
            unsigned int k, unsigned int l,
            unsigned int maximumBlockSize =
            std::numeric_limits<unsigned int>::max()) :
        bemcluster<T>(dofs, op_perm, k, l),
        m_maximumBlockSize(maximumBlockSize) {
    }

    virtual ExtendedBemCluster* clone(unsigned int* op_perm,
                                      unsigned int beg,
                                      unsigned int end) const {
        return new ExtendedBemCluster(cluster_pca<T>::dofs, op_perm, beg, end,
                                      m_maximumBlockSize);
    }

    virtual bool isadm(double eta2, cluster* cl, bl_info& info) {
        if (this->size() > m_maximumBlockSize ||
                cl->size() > m_maximumBlockSize)
            return (info.is_adm = false);
        else
            return bemcluster<T>::isadm(eta2, cl, info);
    }

    unsigned int maximumBlockSize() const {
        return m_maximumBlockSize;
    }

private:
    unsigned int m_maximumBlockSize;
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
boost::shared_array<mblock<typename AhmedTypeTraits<ValueType>::Type>*>
allocateAhmedMblockArray(
        const blcluster* cluster)
{
    typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;
    AhmedMblock** blocks = 0;
    const size_t blockCount = cluster->nleaves();
    allocmbls(blockCount, blocks);
    return boost::shared_array<AhmedMblock*>(
                blocks, AhmedMblockArrayDeleter(blockCount));
}


} // namespace Bempp

#endif
