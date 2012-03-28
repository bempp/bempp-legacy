#ifndef fiber_opencl_options_hpp
#define fiber_opencl_options_hpp

namespace Fiber
{

struct OpenClOptions
{
    OpenClOptions ()
    {
        useOpenCl = true;
    }

    bool useOpenCl;
    // add more as required
};

}

#endif
