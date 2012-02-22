#ifndef fiber_piecewise_linear_continuous_scalar_basis_hpp
#define fiber_piecewise_linear_continuous_scalar_basis_hpp

#include "basis.hpp"

#include "basis_data.hpp"
#include "dune_basis_helper.hpp"

#include <dune/localfunctions/lagrange/p1/p1localbasis.hh>
#include <dune/localfunctions/lagrange/q1/q1localbasis.hh>

namespace Fiber
{

template <int elementVertexCount, typename ValueType>
struct PiecewiseLinearContinuousScalarBasisTraits
{
};

// Line
template <typename ValueType>
struct PiecewiseLinearContinuousScalarBasisTraits<2, ValueType>
{
private:
    typedef ValueType CoordinateType;
public:
    typedef Dune::Q1LocalBasis<CoordinateType, ValueType, 1> DuneBasis;
};

// Triangle
template <typename ValueType>
struct PiecewiseLinearContinuousScalarBasisTraits<3, ValueType>
{
private:
    typedef ValueType CoordinateType;
public:
    typedef Dune::P1LocalBasis<CoordinateType, ValueType, 2> DuneBasis;
};

// Quadrilateral
template <typename ValueType>
struct PiecewiseLinearContinuousScalarBasisTraits<4, ValueType>
{
private:
    typedef ValueType CoordinateType;
public:
    typedef Dune::Q1LocalBasis<CoordinateType, ValueType, 2> DuneBasis;
};

template <int elementVertexCount, typename ValueType>
class PiecewiseLinearContinuousScalarBasis : public Basis<ValueType>
{
private:
    typedef typename PiecewiseLinearContinuousScalarBasisTraits
    <elementVertexCount, ValueType>::DuneBasis DuneBasis;

public:
    virtual int size() const {
        DuneBasis basis;
        return basis.size();
    }

    virtual int order() const {
        return 1;
    }

    virtual void evaluate(int what,
                          const arma::Mat<ValueType>& points,
                          LocalDofIndex localDofIndex,
                          BasisData<ValueType>& data) const
    {
        if (localDofIndex != ALL_DOFS &&
                (localDofIndex < 0 || size() <= localDofIndex))
            throw std::invalid_argument("PiecewiseLinearContinuousScalarBasis::"
                                        "evaluate(): Invalid localDofIndex");

        if (what & VALUES)
            evaluateBasisFunctionsWithDune<ValueType, ValueType, DuneBasis>(
            points, localDofIndex, data.values);
        if (what & DERIVATIVES)
        {
            throw std::runtime_error("PiecewiseLinearContinuousScalarBasis::"
                                     "evaluation of derivatives "
                                     "not implemented yet");
        }
    }
};

} // namespace Fiber

#endif
