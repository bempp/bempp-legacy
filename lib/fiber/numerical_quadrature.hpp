#ifndef fiber_numerical_quadrature_hpp
#define fiber_numerical_quadrature_hpp

/** \file Low-level functions filling arrays of quadrature points and weights. */

#include "element_pair_topology.hpp"

#include <armadillo>
#include <boost/tuple/tuple_comparison.hpp>
#include <ostream>

namespace Fiber
{

struct SingleQuadratureDescriptor
{
    int vertexCount;
    int order;

    bool operator<(const SingleQuadratureDescriptor& other) const {
        return std::make_pair(vertexCount, order) <
                std::make_pair(other.vertexCount, other.order);
    }

    bool operator==(const SingleQuadratureDescriptor& other) const {
        return vertexCount == other.vertexCount &&
                order == other.order;
    }

    bool operator!=(const SingleQuadratureDescriptor& other) const {
        return !operator==(other);
    }

    friend std::ostream&
    operator<< (std::ostream& dest, const SingleQuadratureDescriptor& obj)
    {
        dest << obj.vertexCount << " " << obj.order;
        return dest;
    }
};

struct DoubleQuadratureDescriptor
{
    ElementPairTopology topology;
    int testOrder;
    int trialOrder;

    bool operator<(const DoubleQuadratureDescriptor& other) const {
        using boost::tuples::make_tuple;
        return make_tuple(topology, testOrder, trialOrder) <
                make_tuple(other.topology, other.testOrder, other.trialOrder);
    }

    bool operator==(const DoubleQuadratureDescriptor& other) const {
        return topology == other.topology &&
                testOrder == other.testOrder &&
                trialOrder == other.trialOrder;
    }

    bool operator!=(const DoubleQuadratureDescriptor& other) const {
        return !operator==(other);
    }

    friend std::ostream&
    operator<< (std::ostream& dest, const DoubleQuadratureDescriptor& obj)
    {
        dest << obj.topology << " " << obj.testOrder << " " << obj.trialOrder;
        return dest;
    }
};

template <typename ValueType>
void fillSingleQuadraturePointsAndWeights(int elementCornerCount,
                                          int quadratureOrder,
                                          arma::Mat<ValueType>& points,
                                          std::vector<ValueType>& weights);

template <typename ValueType>
void fillDoubleSingularQuadraturePointsAndWeights(
        const DoubleQuadratureDescriptor& desc,
        arma::Mat<ValueType>& testPoints,
        arma::Mat<ValueType>& trialPoints,
        std::vector<ValueType>& weights);

} // namespace Fiber

#endif
