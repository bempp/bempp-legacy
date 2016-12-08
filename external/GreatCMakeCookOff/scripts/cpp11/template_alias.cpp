template<class A, class B> struct AClass {
  typedef A t_first;
  typedef B t_second;
};

template<class B> using Specialized = AClass<int, B>;

int main() { return 0; }
