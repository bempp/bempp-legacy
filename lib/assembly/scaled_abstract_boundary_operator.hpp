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

#ifndef bempp_scaled_abstract_boundary_operator_hpp
#define bempp_scaled_abstract_boundary_operator_hpp

#include "abstract_boundary_operator_superposition_base.hpp"
#include "boundary_operator.hpp"

namespace Bempp
{

/** \ingroup composite_boundary_operators
 *  \brief Scaled abstract boundary operator.
 *
 *  This class represents an abstract boundary operator multiplied by a scalar. */
template <typename BasisFunctionType_, typename ResultType_>
class ScaledAbstractBoundaryOperator :
        public AbstractBoundaryOperatorSuperpositionBase<BasisFunctionType_, ResultType_>
{
    typedef AbstractBoundaryOperatorSuperpositionBase<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc AbstractBoundaryOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc AbstractBoundaryOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc AbstractBoundaryOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc AbstractBoundaryOperator::QuadratureStrategy */
    typedef typename Base::QuadratureStrategy QuadratureStrategy;
    /** \copydoc AbstractBoundaryOperatorSuperpositionBase::LocalAssembler */
    typedef typename Base::LocalAssembler LocalAssembler;

    /** \brief Constructor.
     *
     *  Construct the boundary operator \f$\alpha L\f$, where
     *  \f$\alpha\f$ is the scalar \p multiplier_ and \f$L\f$ is the operator
     *  represented by \p multiplicand_.
     *
     *  By default the symmetry of the weak form of the resulting operator is
     *  determined automatically. It can be set manually via the parameter \p
     *  symmetry, which can be any combination of the flags defined in the
     *  enumeration type Symmetry. */
    ScaledAbstractBoundaryOperator(
            ResultType multiplier_,
            const BoundaryOperator<BasisFunctionType, ResultType>& multiplicand_,
            int symmetry = AUTO_SYMMETRY);

    virtual bool isLocal() const;

    ResultType_ multiplier() const;
    BoundaryOperator<BasisFunctionType_, ResultType_> multiplicand() const;

private:
    /** \cond PRIVATE */
    ResultType m_multiplier;
    BoundaryOperator<BasisFunctionType, ResultType> m_multiplicand;
    /** \endcond */
};

} // namespace Bempp

#endif
