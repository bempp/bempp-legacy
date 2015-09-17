import bempp
import pickle

class TestPickleGlobalOptions(object):

    def test_pickle_global_options(self):

        serialized = pickle.dumps(bempp.global_parameters)
        restored = pickle.loads(serialized)
        assert bempp.global_parameters.assembly.boundary_operator_assembly_type \
                           == restored.assembly.boundary_operator_assembly_type
        assert bempp.global_parameters.assembly.potential_operator_assembly_type \
                           == restored.assembly.potential_operator_assembly_type
        assert bempp.global_parameters.assembly.enable_singular_integral_caching \
                           == restored.assembly.enable_singular_integral_caching

        assert bempp.global_parameters.quadrature.double_singular \
                           == restored.quadrature.double_singular

        assert bempp.global_parameters.quadrature.near.max_rel_dist \
                           == restored.quadrature.near.max_rel_dist
        assert bempp.global_parameters.quadrature.near.single_order \
                           == restored.quadrature.near.single_order
        assert bempp.global_parameters.quadrature.near.double_order \
                           == restored.quadrature.near.double_order

        assert bempp.global_parameters.quadrature.medium.max_rel_dist \
                           == restored.quadrature.medium.max_rel_dist
        assert bempp.global_parameters.quadrature.medium.single_order \
                           == restored.quadrature.medium.single_order
        assert bempp.global_parameters.quadrature.medium.double_order \
                           == restored.quadrature.medium.double_order

        assert bempp.global_parameters.quadrature.far.single_order \
                           == restored.quadrature.far.single_order
        assert bempp.global_parameters.quadrature.far.double_order \
                           == restored.quadrature.far.double_order

        assert bempp.global_parameters.hmat.assembly_mode \
                           == restored.hmat.assembly_mode
        assert bempp.global_parameters.hmat.min_block_size \
                           == restored.hmat.min_block_size
        assert bempp.global_parameters.hmat.max_block_size \
                           == restored.hmat.max_block_size
        assert bempp.global_parameters.hmat.eta \
                           == restored.hmat.eta
        assert bempp.global_parameters.hmat.eps \
                           == restored.hmat.eps
        assert bempp.global_parameters.hmat.max_rank \
                           == restored.hmat.max_rank
        assert bempp.global_parameters.hmat.compression_algorithm \
                           == restored.hmat.compression_algorithm
        assert bempp.global_parameters.hmat.coarsening \
                           == restored.hmat.coarsening
        assert bempp.global_parameters.hmat.coarseningAccuracy \
                           == restored.hmat.coarseningAccuracy
        assert bempp.global_parameters.hmat.mat_vec_parallel_levels \
                           == restored.hmat.mat_vec_parallel_levels
