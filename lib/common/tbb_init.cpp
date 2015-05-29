
#include "tbb_init.hpp"

namespace Bempp {


    TbbInit::TbbInit() {};

    TbbInit::TbbInit(int number_of_processes) :
        m_task_scheduler(number_of_processes){};


}

