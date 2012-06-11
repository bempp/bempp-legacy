#ifndef bempp_mass_matrix_container_initialiser_hpp
#define bempp_mass_matrix_container_initialiser_hpp

#include <memory>

namespace Bempp
{

template <typename ValueType> class MassMatrixContainer;
template <typename BasisFunctionType> class Space;

template <typename BasisFunctionType, typename ResultType>
class MassMatrixContainerInitialiser
{
public:
    MassMatrixContainerInitialiser(const Space<BasisFunctionType>& space) :
        m_space(space) {
    }

    ~MassMatrixContainerInitialiser();

    std::auto_ptr<MassMatrixContainer<ResultType> > operator()() const;

private:
    const Space<BasisFunctionType>& m_space;
};

} // namespace Bempp

#endif
