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
