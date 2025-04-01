#include <fmt/core.h>

#include "sublibrary/sublib_component.h"

int main() {
    fmt::print("Hello, world!\n");
    int out = name::sublibrary::SubLibComponent::add(2, 3);
    fmt::print("2 + 3 = {}\n", out);
    return 0;
}
