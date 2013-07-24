#ifndef bempp_vc_multicomponent_kernel_values_hpp
#define bempp_vc_multicomponent_kernel_values_hpp

#include "fiber/_2d_array.hpp"

#include <Vc/Vc>

namespace Fiber
{

template <typename ValueType>
struct VcKernelValues
{
     VcKernelValues(size_t testPointCount_, size_t trialPointCount_,
                   size_t rowCount_, size_t colCount_) :
        paddedTestPointCount(((testPointCount_ + Vc::Vector<ValueType>::Size - 1) / 
                              Vc::Vector<ValueType>::Size) * Vc::Vector<ValueType>::Size),
        testPointChunkCount(((testPointCount_ + Vc::Vector<ValueType>::Size - 1) / 
                              Vc::Vector<ValueType>::Size)),
        values(rowCount_, colCount_),
        valuesImag(rowCount_, colCount_)
        {
            for (size_t j = 0; j < values.extent(1); ++j)
                for (size_t i = 0; i < values.extent(0); ++i)
                    values(i, j) = 
                        new Vc::Vector<ValueType>[testPointChunkCount * trialPointCount_];
            for (size_t j = 0; j < values.extent(1); ++j)
                for (size_t i = 0; i < values.extent(0); ++i)
                    valuesImag(i, j) = 
                        new Vc::Vector<ValueType>[testPointChunkCount * trialPointCount_];
            // TODO: destructor
        }

    ~VcKernelValues() {
        for (size_t j = 0; j < values.extent(1); ++j)
            for (size_t i = 0; i < values.extent(0); ++i)
                delete values(i, j);
        for (size_t j = 0; j < values.extent(1); ++j)
            for (size_t i = 0; i < values.extent(0); ++i)
                delete valuesImag(i, j);
    }

    size_t paddedTestPointCount, testPointChunkCount;
    // values(i, j)[k]: kernel row #i, kernel column #j, vector #k
    Fiber::_2dArray<Vc::Vector<ValueType>*> values, valuesImag;
};

}

#endif
