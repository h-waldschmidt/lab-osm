#pragma once

#include <string>

namespace labosm::server {
void simpleServer(const std::string& fmi_file);
void advancedServer(const std::string& fmi_file);
}  // namespace labosm::server