from libcpp.string cimport string

cdef extern from "bempp/assembly/symmetry.hpp":
    cdef enum SymmetryMode "Bempp::Symmetry":
        no_symmetry "Bempp::Symmetry::NO_SYMMETRY"
        symmetric "Bempp::Symmetry::SYMMETRIC"
        hermitian "Bempp::Symmetry::HERMITIAN"
        auto_symmetry "Bempp::Symmetry::AUTO_SYMMETRY"

cdef SymmetryMode symmetry_mode(string name)
