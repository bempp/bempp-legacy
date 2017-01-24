#include <type_traits>

template<int N>
  typename std::enable_if<N==2, int>::type testme() { return 0; }

int main() {
  return testme<2>();
};
