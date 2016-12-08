template<class A, class B=unsigned> B castme(A _a) { 
  return B(_a);
}
int main() {
  unsigned a = castme<int>(-2);
  return 0;
}
