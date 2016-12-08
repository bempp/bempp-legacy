struct A {
  int b;
  A(int _b) : b(_b) {};
  A(A const &) = delete;
};

int main() {
  A first(1);
  return 0;
}
