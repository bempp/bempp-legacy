// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_SCALAR_TRAITS_HPP
#define HMAT_SCALAR_TRAITS_HPP

#include <complex>

namespace hmat {

template <typename T> struct ScalarTraits {

  typedef T RealType;
  typedef T ComplexType;

  ScalarTraits() {
    static_assert(
        sizeof(T) == 0,
        "ScalarTraits only implemented for real/complex float/double types.");
  }
};

template <> struct ScalarTraits<float> {
  typedef float RealType;
  typedef std::complex<float> ComplexType;
};

template <> struct ScalarTraits<double> {
  typedef double RealType;
  typedef std::complex<double> ComplexType;
};

template <> struct ScalarTraits<std::complex<float>> {
  typedef float RealType;
  typedef std::complex<float> ComplexType;
};

template <> struct ScalarTraits<std::complex<double>> {
  typedef double RealType;
  typedef std::complex<double> ComplexType;
};
}

#endif
