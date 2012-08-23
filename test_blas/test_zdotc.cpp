#include <complex>
#include <iostream>

typedef std::complex<double> dcomp;

extern "C" {
dcomp zdotc_(const unsigned*,
             const dcomp*, const unsigned*,
             const dcomp*, const unsigned*);
} // extern "C"

int main()
{
  const unsigned int N = 100;
  const unsigned int ONE = 1;
  dcomp a[N];
  dcomp b[N];

  for (unsigned int i = 0; i < N; ++i) {
      a[i] = dcomp(i, i / 10.);
      b[i] = a[i] / 2.;
  }

  dcomp exact_result(0., 0.);
  for (unsigned int i = 0; i < N; ++i)
    exact_result += conj(a[i]) * b[i];

  dcomp result(0., 0.);
  result = zdotc_(&N, a, &ONE, b, &ONE);
  if (std::abs(result - exact_result) < 1e-13) {
      std::cout << "SUCCESS" << std::endl;
      return 0;
  } else {
      std::cout << "FAILED" << std::endl;
      return 1;
  }
}
