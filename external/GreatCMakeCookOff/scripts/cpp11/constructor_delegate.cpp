#include <cmath>
#include <assert.h>

class A {
  public:
    int a;
    double b;

    A() : A(0, 0) {}
    explicit A(int _a) : A(_a, 0) {}
    explicit A(double _b) : A(0, _b) {}
    A(int _a, double _b) : a(_a), b(_b) {}
};

int main() {
  A instance(1, 1.5e0);
  assert(instance.a == 1);
  assert(std::abs(instance.b - 1.5) < 1e-8);

  A instance1(1);
  assert(instance1.a == 1);
  assert(std::abs(instance1.b) < 1e-8);

  A instance2(1.5);
  assert(instance2.a == 0);
  assert(std::abs(instance2.b - 1.5) < 1e-8);

  A instance3;
  assert(instance3.a == 0);
  assert(std::abs(instance3.b) < 1e-8);
  return 0;
}
