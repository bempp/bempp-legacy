#include <complex>
#include <iostream>


extern "C" {

void daxpy_(const int* n, const double* a, const double* x,
	    const int* incx, double* y, const int* incy);

} // extern "C"

int main()
{
  const int N = 10;
  const double alpha=2.;
  const int ONE = 1;
  double x[N];
  double y[N];
  double res[N];

  for (int i = 0; i < N; ++i) {
      x[i] = 1.0*i;
      y[i] = -2.0*i;
  }

  daxpy_(&N,&alpha,x,&ONE,y,&ONE);

  for (int i = 0;i < N; ++i) {
	  if (std::abs(y[i])>1E-15) return 1;
  }
  std::cout << "Blas run successful" << std::endl;
  return 0;
}
