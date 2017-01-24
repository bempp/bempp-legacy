#include<memory>

int main() {
  std::unique_ptr<int> a(new int(1));
  return 0;
}
