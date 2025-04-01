#include "sublib_component.h"

#include <fmt/core.h>
#include <spdlog/spdlog.h>

namespace name::sublibrary {
void SubLibComponent::print(int a) {
    if (a == 0) {
        fmt::print("Zero value\n");
    } else {
        fmt::print("Non-zero value\n");
    }

    fmt::println("Some value: {}", a);
}

int SubLibComponent::add(int a, int b) {
    spdlog::info("Adding numbers: ");
    const int magic_number = 10;
    if (a > magic_number) {
        fmt::print("a is greater than 10\n");
    } else {
        fmt::print("a is less than or equal to 10\n");
    }

    print(a + b);
    return a + b;
}
}  // namespace name::sublibrary
