#include <iterator>
#include <vector>
#include <exception>

int main() {
  std::vector<int> vector(2, 1);
  auto i_first = std::begin(vector);
  auto i_end = std::end(vector);
  if(i_first + 2 != i_end) return 1;
  if(*i_first != 1) return 2;
  return 0;
}
