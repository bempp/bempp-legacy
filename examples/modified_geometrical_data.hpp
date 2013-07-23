#ifndef bempp_modified_geometrical_data_hpp
#define bempp_modified_geometrical_data_hpp

#include <stdlib.h>

namespace Fiber
{

const int MAX_WORLD_DIM = 3;

template <typename CoordinateType>
struct ModifiedGeometricalData
{
    ModifiedGeometricalData(size_t pointCount_) :
        pointCount(pointCount_) 
        {
            for (int i = 0; i < MAX_WORLD_DIM; ++i) {
                // globals[i] = (CoordinateType*)_mm_malloc(pointCount * sizeof(CoordinateType), 16);
                void* tmp;
                posix_memalign(&tmp, 16, pointCount * sizeof(CoordinateType));

                globals[i] = (CoordinateType*) tmp;
            }
        }

    ~ModifiedGeometricalData()
        {
            for (int i = 0; i < MAX_WORLD_DIM; ++i)
                // _mm_free(globals[i]);
                free(globals[i]);
        }

    size_t pointCount;
    CoordinateType* __restrict __attribute__((aligned)) globals[MAX_WORLD_DIM];
};

}

#endif
