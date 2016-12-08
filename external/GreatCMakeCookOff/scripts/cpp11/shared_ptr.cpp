#include<memory>

int main() {
  std::shared_ptr<int> a(new int(1));
  return 0;
}
