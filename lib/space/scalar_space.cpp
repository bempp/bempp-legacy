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

#include "scalar_space.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/default_collection_of_shapeset_transformations.hpp"

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType>
struct ScalarSpace<BasisFunctionType>::Impl
{
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
    TransformationFunctor;

    Impl() : transformations(TransformationFunctor())
    {}

    Fiber::DefaultCollectionOfShapesetTransformations<TransformationFunctor>
    transformations;
};
/** \endcond */

template <typename BasisFunctionType>
ScalarSpace<BasisFunctionType>::ScalarSpace(const shared_ptr<const Grid>& grid) :
    Base(grid), m_impl(new Impl)
{
}

template <typename BasisFunctionType>
ScalarSpace<BasisFunctionType>::ScalarSpace(
        const ScalarSpace& other) :
    Base(other), m_impl(new Impl(*other.m_impl))
{
}

template <typename BasisFunctionType>
ScalarSpace<BasisFunctionType>::~ScalarSpace()
{
}

template <typename BasisFunctionType>
ScalarSpace<BasisFunctionType>&
ScalarSpace<BasisFunctionType>::
operator=(const ScalarSpace& rhs)
{
    if (this != &rhs) {
        Base::operator=(rhs);
        m_impl.reset(new Impl(*rhs.m_impl));
    }
    return *this;
}

template <typename BasisFunctionType>
const typename ScalarSpace<BasisFunctionType>::CollectionOfShapesetTransformations&
ScalarSpace<BasisFunctionType>::basisFunctionValue() const
{
    return m_impl->transformations;
}

template <typename BasisFunctionType>
void ScalarSpace<BasisFunctionType>::getGlobalDofInterpolationDirections(
        arma::Mat<CoordinateType>& directions) const {
    directions.set_size(1 /* scalar space */, this->globalDofCount());
    directions.fill(1);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(ScalarSpace);

} // namespace Bempp
