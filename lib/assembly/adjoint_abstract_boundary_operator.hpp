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

#ifndef bempp_adjoint_abstract_boundary_operator_hpp
#define bempp_adjoint_abstract_boundary_operator_hpp

#include "abstract_boundary_operator.hpp"
#include "boundary_operator.hpp"

namespace Bempp
{

/** \ingroup composite_boundary_operators
 *  \brief Adjoint abstract boundary operator.
 *
 *  This class represents the adjoint of an abstract boundary operator.
 *
 *  Currently the adjoint operator can be constructed only for
 *  operators acting on real-valued basis functions. Let \f$A\f$ be
 *  such an operator, and let \f$\mathsf{A}\f$ be its discrete weak
 *  form. Then the discrete weak form of the adjoint of \f$A\f$ is the
 *  transpose of \f$\mathsf{A}\f$.
 *
 *  Note that this definition is not the one most usual in mathematics
 *  (where the conjugate transpose would be used). However, it
 *  corresponds to the relation between the double-layer potential
 *  boundary operator and adjoint double-layer potential boundary
 *  operator for Laplace and Helmholtz equation as defined in standard
 *  texts. Hence, for example, the discrete weak forms of the
 *  operators B and C in the code snippet below are numerically
 *  identical:
 *  \code
    BoundaryOperator<BFT, RT> A =
        helmholtz3dDoubleLayerBoundaryOperator<BFT>(
            context, domain, range, dualToRange, waveNumber);
    BoundaryOperator<BFT, RT> B =
        helmholtz3dAdjointDoubleLayerBoundaryOperator<BFT>(
            context, domain, range, dualToRange, waveNumber);
    BoundaryOperator<BFT, RT> C = adjoint(A);
    \endcode
 *  \see adjoint().
 */
template <typename BasisFunctionType_, typename ResultType_>
class AdjointAbstractBoundaryOperator :
        public AbstractBoundaryOperator<BasisFunctionType_, ResultType_>
{
    typedef AbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc AbstractBoundaryOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc AbstractBoundaryOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc AbstractBoundaryOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc AbstractBoundaryOperator::QuadratureStrategy */
    typedef typename Base::QuadratureStrategy QuadratureStrategy;

    /** \brief Constructor.
     *
     *  Construct the boundary operator \f$\alpha L\f$, where
     *  \f$\alpha\f$ is the scalar \p weight and \f$L\f$ is the operator
     *  represented by \p boundaryOp.
     *
     *  By default the symmetry of the weak form of the resulting operator is
     *  determined automatically. It can be set manually via the parameter \p
     *  symmetry, which can be any combination of the flags defined in the
     *  enumeration type Symmetry. */
    AdjointAbstractBoundaryOperator(
            const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
            int symmetry = AUTO_SYMMETRY);

    virtual bool isLocal() const;

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

private:
    BoundaryOperator<BasisFunctionType, ResultType> m_operator;
};

} // namespace Bempp

#endif
