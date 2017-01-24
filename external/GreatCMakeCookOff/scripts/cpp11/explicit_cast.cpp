class A {
  public:
   A(int a = 5) : a_(a) {};
   explicit operator bool() const { return a_ == 2; }
  protected:
   int a_;
};

int main () {
  A a(6);
  A const b(2);
  return (static_cast<bool>(a) == false and static_cast<bool>(b) == true) ? 0: 1;  
}
