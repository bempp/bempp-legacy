
#include "tbb_init.hpp"
#include "py_init.hpp"
#include "eigen_support.hpp"

namespace Bempp {

PyInit PyInit::m_singleton;
EigenInit EigenInit::m_singleton;
TbbInit TbbInit::m_singleton;
}
