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

#ifndef bempp_operator_aca_assembly_helper_hpp
#define bempp_operator_aca_assembly_helper_hpp

//UNFINISHED CLEANUP
//template <typename ValueType>
//class OperatorAcaAssemblyHelper
//{
//public:
//    OperatorAcaAssemblyHelper(const LinearOperator<ValueType>& op,
//              const arma::Mat<ctype>& testPoints,
//              const Space& trialSpace,
//              const arma::Col<unsigned int>& p2oPoints,
//              const arma::Col<unsigned int>& p2oDofs) :
//        m_operator(op), m_grid(grid), m_testPoints(testPoints),
//        m_trialSpace(trialSpace), m_p2oPoints(p2oPoints), m_p2oDofs(p2oDofs)
//    {}

//    //      stores the entries of the block defined
//    //      by b1, n1, b2, n2 (in permuted ordering) in data
//    void cmpbl(unsigned b1, unsigned n1, unsigned b2, unsigned n2,
//               ValueType* data)
//    {
//        arma::Mat<ValueType> dataWrapper(data, n1, n2, false /*copy_aux_mem*/,
//                                         true /*strict*/);

//        arma::Mat<ctype> points(m_testPoints.n_rows, n1);
//        for (int i = 0; i < n1; ++i)
//            points.col(i) = m_testPoints.col(m_p2oPoints(b1 + i));
//        arma::Col<unsigned int> globalDofs(n2);
//        for (int i = 0; i < n2; ++i)
//            globalDofs(i) = m_p2oDofs(b2 + i);

//        // To be decided: who will store the EntityPointers
//        // pointed by LocalDofs?
//        std::vector<LocalDof> localDofs;
//        m_trialSpace.localDofs(globalDofs, localDofs);

//        // Split localDofs into elements: construct
//        //   map<EntityPointer<0>*,
//        //       set<pair<localDofIndex, blockColumn> > >.
//        // (set sorted after localDofIndex)
//        // Then convert to three containers:
//        //   vector<EntityPointer<0>*> elements;
//        //   vector<vector<int> > localDofs;
//        //   vector<vector<int> > blockColumns.
//        //

//        arma::Cube<ValueType> localResult;
//        m_operator.evaluateOperator(m_testPoints, elements, localDofs,
//                                    m_trialSpace, localResult);

//        dataWrapper.zeros();
//        // Iterate over slices of localResult, adding them to appropriate
//        // places in dataWrapper. (Note: it might be difficult for
//        // vector-valued operators unless Ahmed guarantees that all components
//        // will stick together.)

//        // Better to allow only two situations:
//        //   one local dof on one element
//        //   or all local dofs on all elements

//        // For weak form:
//        // list of test elements, list of trial elements (all dofs)
//        // list of test elements, one dof of one trial element
//        // one dof of one test element, list of trial elements
//    }

//    /** Expected size of the entries in this block. */
//    ValueType scale(unsigned b1, unsigned n1, unsigned b2, unsigned n2);

//private:
//    const LinearOperator<ValueType>& m_operator;
//    const Grid& m_grid;
//    const arma::Mat<ctype>& m_testPoints;
//    const Space& m_trialSpace;
//    const arma::Col<unsigned int>& m_p2oPoints;
//    const arma::Col<unsigned int>& m_p2oDofs;
//};

#endif
