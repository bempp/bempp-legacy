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


#ifndef bempp_index_permutation_hpp
#define bempp_index_permutation_hpp

#include <armadillo>

namespace Bempp
{

class IndexPermutation
{
public:
    IndexPermutation(const std::vector<unsigned int>& permutedIndices) :
        m_permutedIndices(permutedIndices) {
    }

    /** \brief Convert a vector from original to permuted ordering. */
    template <typename ValueType>
    void permuteVector(const arma::Col<ValueType>& original,
                       arma::Col<ValueType>& permuted) const
    {
        const int dim = original.n_elem;
        permuted.set_size(dim);
        for (int i = 0; i < dim; ++i)
            permuted(m_permutedIndices[i]) = original(i);
    }

    /** \brief Convert a vector from permuted to original ordering. */
    template <typename ValueType>
    void unpermuteVector(const arma::Col<ValueType>& permuted,
                         arma::Col<ValueType>& original) const
    {
        const int dim = permuted.n_elem;
        original.set_size(dim);
        for (int i = 0; i < dim; ++i)
            original(i) = permuted(m_permutedIndices[i]);
    }

    /** \brief Permute index. */
    unsigned int permuted(unsigned int index) const {
        return m_permutedIndices[index];
    }

private:
    const std::vector<unsigned int> m_permutedIndices;
};

} // namespace Bempp

#endif
