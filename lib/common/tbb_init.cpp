#include <boost/lexical_cast.hpp>
#include "tbb_init.hpp"
#include <cstdlib>

namespace Bempp {

TbbInit::TbbInit() {

  unsigned int num_threads = 0;
  char *p_num_threads = std::getenv("BEMPP_NUM_THREADS");
  if (p_num_threads != nullptr) {
    try {
      num_threads = boost::lexical_cast<unsigned int>(p_num_threads);
    } catch (const boost::bad_lexical_cast &) {
    }
  }
  if (num_threads < 1) {
    m_task_scheduler.reset(new tbb::task_scheduler_init());

  } else {
    m_task_scheduler.reset(new tbb::task_scheduler_init(num_threads));
  }
}
}
