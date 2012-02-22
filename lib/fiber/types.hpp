#ifndef fiber_types_hpp
#define fiber_types_hpp

namespace Fiber
{

enum CallVariant
{
    TEST_TRIAL = 0,
    TRIAL_TEST = 1
};

typedef int LocalDofIndex;
const LocalDofIndex ALL_DOFS = -1;

}

#endif
