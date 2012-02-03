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

#ifndef bempp_space_hpp
#define bempp_space_hpp

#include "../common/types.hpp"

namespace Bempp
{

template <int codim> class Entity;

template <typename ValueType>
class Space
{
public:
    // A grid reference is necessary because e.g. when setting the element
    // variant it is necessary to check whether the element is triangular
    // or quadrilateral. Also, requests for element refinement should probably
    // be made via Space rather than via Grid.
    Space(Grid& grid) : m_grid(grid)
    {}

    /** @name Attributes
    @{ */

    virtual int domainDimension() const = 0;
    /** Number of components (e.g. H1 space -> 1, H(curl) space on a 2D surface -> 2) */
    virtual int codomainDimension() const = 0;
    virtual int basisFunctionCount(ElementVariant elementVariant) const = 0;

    const Grid& grid() const { return m_grid; }

    /** @}
        @name Function evaluation
        @{ */
    virtual void evaluateBasisFunctions(
            // interpreted differently for different spaces.
            // e.g. it could be a bitfield
            // {4 /* quad */, 3 /* polynomial order in x dir. */,
            //  2 /* polynomial order in y dir. */}
            ElementVariant elementVariant,
            // rows: local coordinates, cols: points
            const arma::Mat<ctype>& local,
            /*
            // which functions to evaluate?
            std::vector<int> functionIds,
            */
            // rows: components, cols: points, slices: basis function
            arma::Cube<ValueType>& result) const = 0;
    virtual void evaluateBasisFunctionDerivative(
            // interpreted differently for different spaces.
            // e.g. it could be a bitfield
            // {4 /* quad */, 3 /* polynomial order in x dir. */,
            //  2 /* polynomial order in y dir. */}
            ElementVariant elementVariant,
            // rows: points, cols: individual (x, y) coordinates
            const arma::Mat<ctype>& local,
            /*
            // which functions to evaluate?
            std::vector<int> functionIds,
            */
            // local coordinate along which the derivative is taken
            int direction,
            // rows: components, cols: points, slices: basis function derivatives
            arma::Cube<ValueType>& result) const = 0;

    virtual bool shapeFunctionsDependOnJacobianMatrix() const = 0;

    // result[codomainDimension,basisFunctionCount(),pointCount]
    virtual void evaluateShapeFunctions(
            const EntityPointer<0>* element,
            const arma::Mat<ctype>& local,
            arma::Cube<ValueType>& result) const = 0;
    virtual void evaluateShapeFunctionsInternal(
            const arma::Cube<ValueType>& basisFunctionValues,
            const arma::Cube<ValueType>& jacobianInverseTransposed,
            arma::Cube<ValueType>& result) const = 0;
    // variants for only selected local DOFs also possible.

    virtual void evaluateShapeFunctionSurfaceCurls(
            const EntityPointer<0>* element,
            const arma::Mat<ctype>& local,
            arma::Cube<ValueType>& result) const = 0;
    virtual void evaluateShapeFunctionSurfaceCurlsInternal(
            const arma::Cube<ValueType>& basisFunctionValues,
            const arma::Cube<ValueType>& jacobianInverseTransposed,
            arma::Cube<ValueType>& result) const = 0;

    virtual std::string openClCodeToEvaluateBasisFunctions() const = 0;
    virtual std::string openClCodeToEvaluateShapeFunctions() const = 0;
    virtual std::string openClCodeToEvaluateShapeFunctionSurfaceCurls() const = 0;

    /** @}
        @name Element order management
        @{ */
    virtual void setElementVariant(const Entity<0>& element, ElementVariant variant) = 0;
    virtual ElementVariant elementVariant(const Entity<0>& element) const = 0;
    // additional functions for e.g. increasing polynomial order of all elements
    // ...

    /** @}
        @name DOF management
        @{ */

    void assignDofs();
    virtual bool dofsAssigned() const = 0; // returns a flag that is set to true via assignDofs()

    virtual int globalDofCount() const = 0;
    virtual void globalDofs(const Entity<0>& element,
                            std::vector<GlobalDofIndex>& dofs) const = 0;

    virtual void global2localDofs(
            const std::vector<GlobalDofIndex>& globalDofs,
            std::vector<std::vector<LocalDof> >& localDofs) const = 0;

    // This function is used only by the ACA assembler.
    // For the moment, Point will always be 3D, independently from the
    // actual dimension of the space. Once Ahmed's bemcluster is made dimension-
    // independent, we may come up with a more elegant solution.
    virtual void globalDofPositions(arma::Col<Point3D>& positions) const = 0;
    /** @} */

private:
    virtual void assignVertexDofs() = 0;
    virtual void assignEdgeDofs() = 0;
    virtual void assignBubbleDofs() = 0;

protected:
    Grid& m_grid;
};

} //namespace Bempp

#endif
