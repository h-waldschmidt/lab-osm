#include "sublib_component.h"

#include <gtest/gtest.h>

namespace name::sublibrary {
TEST(SubLibComponentTest, Add) {
    int out = name::sublibrary::SubLibComponent::add(2, 3);
    EXPECT_EQ(out, 5);
}
}  // namespace name::sublibrary
