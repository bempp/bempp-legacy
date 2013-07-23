#ifndef bempp_modified_kernel_values_hpp
#define bempp_modified_kernel_values_hpp

#include <stdlib.h>
#include <iostream>

namespace Fiber
{

template <typename ValueType>
struct ModifiedKernelValues
{
    ModifiedKernelValues(int rowCount_, int colCount_, size_t testPointCount_, size_t trialPointCount_)
        {
            rowCount = rowCount_;
            colCount = colCount_;
            values = new ValueType* __restrict *[rowCount];
            for (int i = 0; i < rowCount; ++i) {
                values[i] = new ValueType*[colCount];
                for (int j = 0; j < colCount; ++j) {
                    // values[i][j] = (ValueType*)_mm_malloc(testPointCount_ * trialPointCount_ * sizeof(ValueType), 16);
                    void* tmp;
                    posix_memalign(&tmp, 16, testPointCount_ * trialPointCount_ * sizeof(ValueType));
                    values[i][j] = (ValueType*) tmp;
                }
}
        }

    ~ModifiedKernelValues()
        {
            for (int i = 0; i < rowCount; ++i)
                for (int j = 0; j < colCount; ++j)
                    free(values[i][j]);
        }

    int rowCount, colCount;
    ValueType* __restrict ** values;
};

}

#endif
