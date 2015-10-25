from libcpp.string cimport string

cdef extern from "bempp/assembly/symmetry.hpp":
    cdef enum SymmetryMode "Bempp::Symmetry":
        no_symmetry "Bempp::Symmetry::NO_SYMMETRY"
        symmetric "Bempp::Symmetry::SYMMETRIC"
        hermitian "Bempp::Symmetry::HERMITIAN"
        auto_symmetry "Bempp::Symmetry::AUTO_SYMMETRY"

cdef extern from "bempp/assembly/transposition_mode.hpp":
    cdef enum TranspositionMode "Bempp::TranspositionMode":
        no_transpose "Bempp::TranspositionMode::NO_TRANSPOSE"
        conjugate "Bempp::TranspositionMode::CONJUGATE"
        transpose "Bempp::TranspositionMode::TRANSPOSE"
        conjugate_transpose "Bempp::TranspositionMode::CONJUGATE_TRANSPOSE"

cdef extern from "bempp/hmat/common.hpp":
    cdef enum HMatBlockType "hmat::DataBlockType":
        dense "hmat::DataBlockType::DENSE"
        low_rank_ab "hmat::DataBlockType::LOW_RANK_AB"


cdef SymmetryMode symmetry_mode(string name)
cdef TranspositionMode transposition_mode(string name)
cdef HMatBlockType hmat_block_type(string name)
cdef TranspositionMode compute_transpose_mode(
        TranspositionMode current_mode,
        TranspositionMode input_mode)
