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

#include "../common/common.hpp"

// to ensure there are no inconsistent forward declarations
#include "ahmed_aux_fwd.hpp"

#include "ahmed_leaf_cluster_array.hpp"
#include "ahmed_mblock_array_deleter.hpp"
#include "../common/types.hpp"
#include "../common/bounding_box.hpp"

#include "../common/boost_scoped_array_fwd.hpp"
#include "../common/boost_shared_array_fwd.hpp"
#include <complex>
#include <iostream>
#include <memory>
#include <stdexcept>

// Ahmed's include files

#ifdef __INTEL_COMPILER
#pragma warning(disable:381)
#endif

#include "ahmed_complex.hpp"
#include "extended_bem_cluster.hpp"

#include <bbxbemblcluster.h>
#include <bllist.h>
// #include <matgen_sqntl.h>
// #define _OPENMP
// #include <matgen_omp.h>
// #undef _OPENMP

bool multaHvec_omp(double d, blcluster* bl, mblock<double>** A, double* x,
           double* y);

#include <matrix.h>
#undef SIGN
#undef SQR
#undef MIN
#undef MAX

template <typename T>
T MIN(T x, T y)
{
    return std::min(x, y);
}

template <typename T>
T MAX(T x, T y)
{
    return std::max(x, y);
}

#ifdef __INTEL_COMPILER
#pragma warning(default:381)
#endif


namespace Bempp
{

/** \brief An Ahmed-compatible degree-of-freedom type. */
template <typename CoordinateType>
// struct AhmedDofWrapper : public Point3D<CoordinateType>
struct AhmedDofWrapper : public BoundingBox<CoordinateType>
{
    // Methods required by the AHMED clustering code

    CoordinateType getcenter(int i) const {
        if (i == 0) return this->reference.x;
        else if (i == 1) return this->reference.y;
        else return this->reference.z;
    }

    CoordinateType getlbound(int i) const {
        if (i == 0) return this->lbound.x;
        else if (i == 1) return this->lbound.y;
        else return this->lbound.z;
    }

    CoordinateType getubound(int i) const {
        if (i == 0) return this->ubound.x;
        else if (i == 1) return this->ubound.y;
        else return this->ubound.z;
    }
};

template <typename ValueType>
boost::shared_array<mblock<typename AhmedTypeTraits<ValueType>::Type>*>
allocateAhmedMblockArray(size_t blockCount)
{
    typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;
    AhmedMblock** blocks = 0;
    allocmbls(blockCount, blocks);
    return boost::shared_array<AhmedMblock*>(
                blocks, AhmedMblockArrayDeleter(blockCount));
}

template <typename ValueType>
boost::shared_array<mblock<typename AhmedTypeTraits<ValueType>::Type>*>
allocateAhmedMblockArray(const blcluster* cluster)
{
    return allocateAhmedMblockArray<ValueType>(cluster->nleaves());
}

} // namespace Bempp

#endif
