#include "ahmed_mblock_array_deleter.hpp"
#include "ahmed_aux.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename ValueType>
void AhmedMblockArrayDeleter::operator () (mblock<ValueType>** blocks) const
{
    freembls(m_arraySize, blocks);
}

#define INSTANTIATE_OPERATOR_CALL(RESULT) \
    template void AhmedMblockArrayDeleter::operator () ( \
        mblock<AhmedTypeTraits<RESULT>::Type>** blocks) const

FIBER_ITERATE_OVER_VALUE_TYPES(INSTANTIATE_OPERATOR_CALL);

} // namespace Bempp
