

def compression_statistics(hmatrix):

    leaf_nodes = hmatrix.block_cluster_tree.leaf_nodes

    number_of_dense_blocks = 0
    number_of_low_rank_blocks = 0

    overall_memory_in_kb = 0

    dtype_size = 0

    if hmatrix.dtype=='float64':
        dtype_size = 8*1.0/1024
    else:
        dtype_size = 16*1.0/1024

    dense_memory_in_kb = dtype_size*hmatrix.shape[0]*hmatrix.shape[1]

    for leaf in leaf_nodes:
        data = hmatrix.data(leaf)
        if data.block_type=='dense':
            number_of_dense_blocks += 1
        else:
            number_of_low_rank_blocks += 1
        overall_memory_in_kb += data.mem_size

    total_number_of_blocks = number_of_dense_blocks+number_of_low_rank_blocks

    return {'dense_blocks':number_of_dense_blocks,
            'low_rank_blocks':number_of_low_rank_blocks,
            'ratio_dense_blocks':(1.0*number_of_dense_blocks)/total_number_of_blocks,
            'ratio_low_rank_blocks':(1.0*number_of_low_rank_blocks)/total_number_of_blocks,
            'mem_size_in_mb':(1.0*overall_memory_in_kb)/1024,
            'compression_rate':(1.0*overall_memory_in_kb)/dense_memory_in_kb}




    

