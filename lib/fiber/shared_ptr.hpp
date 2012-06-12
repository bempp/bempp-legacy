#ifndef fiber_shared_ptr_hpp
#define fiber_shared_ptr_hpp

#include "../common/common.hpp"


#ifdef __INTEL_COMPILER
#pragma warning(disable:279 858)
#endif

#include <boost/shared_ptr.hpp>

#ifdef __INTEL_COMPILER
#pragma warning(default:858)
#endif

namespace Fiber
{

using boost::shared_ptr;

struct null_deleter
{
    void operator() (const void*) const {
    }
};

/** \brief Create a shared pointer from a reference to an object allocated on stack.

The object will not be deleted when the shared pointer goes out of scope. */
template <typename T>
inline shared_ptr<T> make_shared_from_ref(T& t)
{
    return shared_ptr<T>(&t, null_deleter());
}

} // namespace Fiber

#endif // fiber_shared_ptr_hpp
