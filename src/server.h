#pragma once

#include <string>

namespace labosm::server {

/**
 * @brief Starts a simple server that only supports Dijkstra for routing.
 * @param fmi_file The path to the fmi graph file.
 */
void simpleServer(const std::string& fmi_file);

/**
 * @brief Starts an advanced server that supports CH and Hub Labeling for routing.
 * @param fmi_file The path to the fmi graph file.
 */
void advancedServer(const std::string& fmi_file);
}  // namespace labosm::server