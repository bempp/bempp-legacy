#include <random>
#include <sstream>

int main()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 6);
    std::ostringstream sstr;
    for(int n=0; n<10; ++n) sstr << dis(gen) << ' ';
    return 0;
}

