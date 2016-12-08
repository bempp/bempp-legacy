#include <cmath>
#include <assert.h>

class A {
  public:
    int a;
    explicit A(int _a) : A(_a) {}
};
class B : public A {
    public:
    using A::A;
};

int main() {
  B instance(2);
  return 0;
}
