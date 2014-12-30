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

cdef extern from "bempp/assembly/grid_function.hpp":
    cdef enum ConstructionMode "Bempp::ConstructionMode":
        approximate "Bempp::ConstructionMode::APPROXIMATE"
        interpolate "Bempp::ConstructionMode::INTERPOLATE"

cdef SymmetryMode symmetry_mode(string name)
cdef TranspositionMode transposition_mode(string name)
cdef ConstructionMode construction_mode(string name)
