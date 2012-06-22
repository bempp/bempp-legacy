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

#ifndef bempp_mass_matrix_container_initialiser_hpp
#define bempp_mass_matrix_container_initialiser_hpp

#include "../common/common.hpp"
#include "config_trilinos.hpp"

#include <memory>

namespace Bempp
{

template <typename ValueType> class MassMatrixContainer;
template <typename BasisFunctionType> class Space;

template <typename BasisFunctionType, typename ResultType>
class MassMatrixContainerInitialiser
{
public:
    MassMatrixContainerInitialiser(const Space<BasisFunctionType>& space,
                                   const Space<BasisFunctionType>& dualSpace,
                                   bool forceDense = false);

    std::auto_ptr<MassMatrixContainer<ResultType> > operator()() const;

private:
    std::auto_ptr<MassMatrixContainer<ResultType> >
    assembleOperatorsInDenseMode() const;
#ifdef WITH_TRILINOS
    std::auto_ptr<MassMatrixContainer<ResultType> >
    assembleOperatorsInSparseMode() const;
#endif

private:
    const Space<BasisFunctionType>& m_space;
    const Space<BasisFunctionType>& m_dualSpace;
    bool m_forceDense;
};

} // namespace Bempp

#endif
