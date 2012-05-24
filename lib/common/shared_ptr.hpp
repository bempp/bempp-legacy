#ifndef bempp_shared_ptr_hpp
#define bempp_shared_ptr_hpp

#include "../fiber/shared_ptr.hpp"

namespace Bempp
{

using Fiber::shared_ptr;
using Fiber::make_shared_from_ref;
using Fiber::null_deleter;

} // namespace Bempp

#endif // bempp_shared_ptr_hpp
