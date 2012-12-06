#include <iostream>
#include <cmath>

extern "C" {
  float sdot_(const unsigned*, const float*, const unsigned*,
              const float*, const unsigned*);
} // extern "C"

int main()
{
  const unsigned int N = 100;
  const unsigned int ONE = 1;
  float a[N];
  float b[N];

  for (unsigned int i = 0; i < N; ++i) {
      a[i] = float(i);
      b[i] = a[i] / 2.;
  }

  float exact_result = 0.;
  for (unsigned int i = 0; i < N; ++i)
    exact_result += a[i] * b[i];

  float result = sdot_(&N, a, &ONE, b, &ONE);
  if (std::abs(result - exact_result) < 1e-6 * fabs(exact_result)) {
      std::cout << "SUCCESS" << std::endl;
      return 0;
  } else {
      std::cout << "FAILED" << std::endl;
      return 1;
  }
}
