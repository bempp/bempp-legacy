#ifndef bempp_ahmed_mblock_array_deleter_hpp
#define bempp_ahmed_mblock_array_deleter_hpp

#include "../common/common.hpp"
#include "ahmed_aux_fwd.hpp"

#include <cstddef>

namespace Bempp
{

class AhmedMblockArrayDeleter
{
public:
    explicit AhmedMblockArrayDeleter(size_t arraySize) :
        m_arraySize(arraySize) {
    }

    template <typename ValueType>
    void operator() (mblock<ValueType>** blocks) const;

private:
    size_t m_arraySize;
};

} // namespace Bempp

#endif
