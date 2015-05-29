#ifndef tbb_init_hpp
#define tbb_init_hpp

#include <tbb/task_scheduler_init.h>

namespace Bempp {

class TbbInit {
    
private:

    TbbInit();
    TbbInit(int number_of_processes);
    TbbInit(const TbbInit & other);
    const TbbInit &operator=(const TbbInit &other);

    static TbbInit m_singleton;

    tbb::task_scheduler_init m_task_scheduler;


};
}

#endif
