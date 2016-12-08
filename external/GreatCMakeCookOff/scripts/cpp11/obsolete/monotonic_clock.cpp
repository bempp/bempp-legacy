#include <chrono>

int main(int argc, char const * argv[]) {
    auto const start = std::chrono::monotonic_clock::now();
    auto const end = std::chrono::monotonic_clock::now();
    // If this works as intended start is always smaller than end.
    // Return 0 i.e. false to indicate success.
    return start > end;
}
