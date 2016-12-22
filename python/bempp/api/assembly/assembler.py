"""Various assemblers for integral operators."""


class LocalOperatorLocalAssembler(object):
    """This assembler evaluates local weak forms for local operators."""

    def __init__(self, impl):
        """Constructor. Should not be called directly."""
        self._impl = impl

    def evaluate_local_weak_forms(self, element_indices):
        """Return local element matrices on the given element indices."""
        return self._impl.evaluate_local_weak_forms(element_indices)

    @property
    def dtype(self):
        """The data type ('float64' or 'complex128') of the local assembler."""
        return self._impl.dtype()


class IntegralOperatorLocalAssembler(object):
    """This assembler evaluates local weak forms for non-local operators."""

    def __init__(self, impl):
        """Constructor. Should not be called directly."""
        self._impl = impl

    def evaluate_local_weak_forms(self, index_pairs):
        """
        Return local element matrices for a list of index pairs.

        Each index is of the form (test_index, trial_index).

        """
        test_indices = [ind_pair[0] for ind_pair in index_pairs]
        trial_indices = [ind_pair[1] for ind_pair in index_pairs]

        return self._impl.evaluate_local_weak_forms(
            test_indices, trial_indices)

    @property
    def dtype(self):
        """The data type ('float64' or 'complex128') of the local assembler."""
        return self._impl.dtype()


def assemble_dense_block(operator, rows, cols, parameters=None):
    """Assemble a dense (sub)-block of an elementary integral operator."""
    from .boundary_operator import ElementaryBoundaryOperator
    from .discrete_boundary_operator import DenseDiscreteBoundaryOperator
    from bempp.core.assembly.assembler import assemble_dense_block_ext

    if not isinstance(operator, ElementaryBoundaryOperator):
        raise TypeError(
            "operator must be of type 'ElementaryBoundaryOperator.")

    if parameters is None:
        parameters = operator.parameters

    return DenseDiscreteBoundaryOperator(
        assemble_dense_block_ext(rows, cols,
                                 operator.domain._impl,
                                 operator.dual_to_range._impl,
                                 operator.local_assembler._impl,
                                 parameters).as_matrix())


def assemble_singular_part(operator):
    """Assemble the singular part of an integral operator."""
    def create_sparse_matrix_from_integration_pairs(
            all_test_trial_function_pairs,
            test_dof_map, trial_dof_map, test_dof_count, trial_dof_count):

        data = np.zeros(len(all_test_trial_function_pairs),
                        dtype=assembler.dtype)
        row_indices = np.zeros(len(all_test_trial_function_pairs),
                               dtype=np.int)
        col_indices = np.zeros(len(all_test_trial_function_pairs),
                               dtype=np.int)

        for index, (test_index, trial_index) in enumerate(
                all_test_trial_function_pairs):

            row_indices[index] = test_index
            col_indices[index] = trial_index

            # Get local dofs and dof_weights
            test_dofs, test_weights = (test_dof_map[0][test_index],
                                       test_dof_map[1][test_index])
            trial_dofs, trial_weights = (trial_dof_map[0][trial_index],
                                         trial_dof_map[1][trial_index])

            for i, test_dof in enumerate(test_dofs):
                for j, trial_dof in enumerate(trial_dofs):
                    element_pair = (test_dof.entity_index,
                                    trial_dof.entity_index)
                    weak_form_data = weak_form_lookup.get(element_pair)
                    data[index] += (weak_form_data[test_dof.dof_index,
                                    trial_dof.dof_index] * np.conj(
                                     test_weights[i]) * np.conj(
                                      trial_weights[j]))

        return SparseDiscreteBoundaryOperator(
            csc_matrix((data, (row_indices, col_indices)),
                       shape=(test_dof_count, trial_dof_count)))

    from scipy.sparse import csc_matrix
    from bempp.api.assembly.discrete_boundary_operator import \
        SparseDiscreteBoundaryOperator
    import numpy as np

    test_space = operator.dual_to_range
    trial_space = operator.domain

    test_dof_count = operator.dual_to_range.global_dof_count
    trial_dof_count = operator.domain.global_dof_count

    # If the test and trial grid are different return a zero matrix

    if operator.domain.grid != operator.dual_to_range.grid:
        return SparseDiscreteBoundaryOperator(
            csc_matrix((test_dof_count, trial_dof_count)))

    # Cache the global to local accesses
    test_dof_map = test_space.global_to_local_dofs(range(test_dof_count))
    trial_dof_map = trial_space.global_to_local_dofs(range(trial_dof_count))

    grid = operator.domain.grid

    # Now get adjacent element pairs

    vertex_to_element_matrix = grid.leaf_view.vertex_to_element_matrix
    element_to_element_matrix = vertex_to_element_matrix.transpose().dot(
        vertex_to_element_matrix)

    nonzero_pairs = element_to_element_matrix.nonzero()
    index_pairs = zip(nonzero_pairs[0], nonzero_pairs[1])

    # Now get all pairs of basis functions who partially
    # overlap via adjacent elements

    all_test_trial_function_pairs = []

    for pair in index_pairs:

        test_element = grid.leaf_view.element_from_index(pair[0])
        trial_element = grid.leaf_view.element_from_index(pair[1])

        global_test_dofs = test_space.get_global_dofs(test_element)
        global_trial_dofs = trial_space.get_global_dofs(trial_element)

        for test_dof_index in global_test_dofs:
            if test_dof_index > -1:
                for trial_dof_index in global_trial_dofs:
                    if trial_dof_index > -1:
                        all_test_trial_function_pairs.append(
                            (test_dof_index, trial_dof_index))

    # Remove duplicates
    all_test_trial_function_pairs = list(set(all_test_trial_function_pairs))

    # Now get all integraton element pairs associated

    all_integration_element_pairs = []

    for function_pair in all_test_trial_function_pairs:

        test_local_dofs = test_dof_map[0][function_pair[0]]
        trial_local_dofs = trial_dof_map[0][function_pair[1]]

        for test_dof in test_local_dofs:
            for trial_dof in trial_local_dofs:
                all_integration_element_pairs.append(
                        (test_dof.entity_index, trial_dof.entity_index))

    # Remove duplicates
    all_integration_element_pairs = list(set(all_integration_element_pairs))

    # Now compute all local dof interactions

    assembler = operator.local_assembler

    weak_forms = assembler.evaluate_local_weak_forms(
        all_integration_element_pairs)

    # Create a dictionary

    weak_form_lookup = dict(zip(all_integration_element_pairs, weak_forms))

    # Now need to create the sparse matrix

    return create_sparse_matrix_from_integration_pairs(
            all_test_trial_function_pairs,
            test_dof_map, trial_dof_map, test_dof_count, trial_dof_count)
