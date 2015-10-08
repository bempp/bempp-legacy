from unittest import TestCase
import bempp
import bempp.api
import pickle


class TestPickleGlobalOptions(TestCase):

    def test_pickle_global_options(self):
        def _iterateorcompare(optionsdict, original, restored):
            for key in optionsdict.keys():
                suboriginal = getattr(original, key)
                subrestored = getattr(restored, key)
                suboptionsdict = optionsdict[key]
                for element in suboptionsdict:
                    if isinstance(element, str):
                        originalval = getattr(suboriginal, element)
                        restoredval = getattr(subrestored, element)
                        assert originalval == restoredval
                        print("testing " + str(originalval))
                    else:
                        _iterateorcompare(element, suboriginal, subrestored)

        serialized = pickle.dumps(bempp.api.global_parameters)
        original = bempp.api.global_parameters
        restored = pickle.loads(serialized)
        optionsdict = {}
        optionsdict['assembly'] = ['boundary_operator_assembly_type',
                                   'potential_operator_assembly_type',
                                   'enable_singular_integral_caching']
        quadraturesuboptions = {}
        quadraturesuboptions['far'] = ['single_order', 'double_order']
        quadraturesuboptions['medium'] = ['max_rel_dist', 'single_order',
                                          'double_order']
        quadraturesuboptions['near'] = quadraturesuboptions['medium']
        optionsdict['quadrature'] = ['double_singular',
                                     quadraturesuboptions]
        optionsdict['hmat'] = ['admissibility', 'coarsening',
                               'coarseningAccuracy', 'compression_algorithm',
                               'eps', 'eta', 'mat_vec_parallel_levels',
                               'max_block_size', 'max_rank']
        _iterateorcompare(optionsdict, original, restored)


if __name__ == "__main__":
    from unittest import main

    main()
