class A {
  public:
    virtual void foo(int &_a) = 0;
    void bar(int &_a)  { _a += 1; }
};

class B : public A {
  public:
    void foo(int &_a) override { _a += 1; }
    void bar(int &_a) override { _a += 2; }
};

int main() {
  A * b = static_cast<A*>(new B);
  int a = 0;
  b->foo(a); 
  if(a != 1) return 1;
  b->bar(a); 
  if(a != 3) return 1;
  delete b;
  return 0;
}
