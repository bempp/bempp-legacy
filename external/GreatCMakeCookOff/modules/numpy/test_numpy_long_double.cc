#include <Python.h>
#include <numpy/npy_common.h>

template<class T> class type;

template<> struct type<npy_longdouble>
{
  typedef npy_longdouble np_type;
  static const int value;
};
const int type<npy_longdouble>::value = 0;
template<> struct type<npy_double>
{
  typedef npy_double np_type;
  static const int value;
};
const int type<npy_double>::value = 0;
int main() {return 0;}
