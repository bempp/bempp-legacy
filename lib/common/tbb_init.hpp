#ifndef tbb_init_hpp
#define tbb_init_hpp

#include <tbb/task_scheduler_init.h>
#include "shared_ptr.hpp"

namespace Bempp {

class TbbInit {

private:
  TbbInit();
  TbbInit(const TbbInit &other);
  const TbbInit &operator=(const TbbInit &other);

  static TbbInit m_singleton;

  shared_ptr<tbb::task_scheduler_init> m_task_scheduler;
};
}

#endif
