#include <Python.h>
#include <numpy/npy_common.h>

template<class T> class type;

template<> struct type<npy_bool>
{
  typedef npy_bool np_type;
  static const int value;
};
const int type<npy_bool>::value = 0;
template<> struct type<npy_ubyte>
{
  typedef npy_ubyte np_type;
  static const int value;
};
const int type<npy_ubyte>::value = 0;
int main() {return 0;}
