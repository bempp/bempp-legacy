#ifndef py_hmat_support_hpp
#define py_hmat_support_hpp

#include "bempp/hmat/common.hpp"
#include "bempp/hmat/hmatrix_low_rank_data.hpp"
#include "bempp/hmat/hmatrix_dense_data.hpp"

namespace hmat {


    template <typename ValueType>
    shared_ptr<const HMatrixLowRankData<ValueType>>
    py_cast_to_low_rank_data(
            const shared_ptr<const HMatrixData<ValueType>>& data)
    {
        return static_pointer_cast<const HMatrixLowRankData<ValueType>>(
                data);
    }

    template <typename ValueType>
    shared_ptr<const HMatrixDenseData<ValueType>>
    py_cast_to_dense_data(
            const shared_ptr<const HMatrixData<ValueType>>& data)
    {
        return static_pointer_cast<const HMatrixDenseData<ValueType>>(
                data);
    }


}


#endif 
