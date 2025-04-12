#pragma once

#include <cstdint>
#include <vector>

namespace labosm {
struct Edge {
    uint32_t m_target;         // Target node ID
    uint32_t m_cost;           // Cost/Distance of the edge // TODO: Rename for consistency
    uint32_t m_child_1_index;  // index to underlying edge if shortcut edge
    uint32_t m_child_2_index;  // index to underlying edge if shortcut edge

    Edge(uint32_t target, uint32_t cost) : m_target(target), m_cost(cost), m_child_1_index(-1), m_child_2_index(-1) {}
    Edge(uint32_t target, uint32_t cost, uint32_t child_1, uint32_t child_2)
        : m_target(target), m_cost(cost), m_child_1_index(child_1), m_child_2_index(child_2) {}

    bool isShortcut() { return (m_child_1_index != -1) && (m_child_2_index != -1); }
};

struct ContractionEdge {
    uint32_t m_target;  // Target node ID
    uint32_t m_cost;    // Cost/Distance of the edge // TODO: Rename for consistency
    uint32_t contracted_node_id;
    uint32_t m_child_1_index;         // index to underlying edge if shortcut edge
    uint32_t m_child_2_index;         // index to underlying edge if shortcut edge
    uint32_t m_num_underlying_edges;  // number of underlying edges (used for MIXED heuristic)

    ContractionEdge(uint32_t target, uint32_t cost)
        : m_target(target), m_cost(cost), contracted_node_id(-1), m_child_1_index(-1), m_child_2_index(-1) {}
    ContractionEdge(uint32_t target, uint32_t cost, uint32_t contracted_node, uint32_t child_1, uint32_t child_2,
                    uint32_t num_underlying_edges)
        : m_target(target),
          m_cost(cost),
          contracted_node_id(contracted_node),
          m_child_1_index(child_1),
          m_child_2_index(child_2),
          m_num_underlying_edges(num_underlying_edges) {}

    bool isShortcut() const { return (m_child_1_index != -1) && (m_child_2_index != -1); }
};
struct ContractionData {
    std::vector<bool> m_outgoing;
    std::vector<uint32_t> m_reset_outgoing;
    std::vector<bool> m_visited;
    std::vector<uint32_t> m_reset_visited;
    std::vector<uint32_t> m_distances;
    std::vector<uint32_t> m_reset_distances;
    std::vector<uint32_t> m_num_contracted_neighbors;
    std::vector<std::pair<uint32_t, ContractionEdge>> m_shortcuts_fwd;
    std::vector<std::pair<uint32_t, ContractionEdge>> m_shortcuts_bwd;

    ContractionData(uint32_t num_nodes)
        : m_outgoing(num_nodes, false),
          m_visited(num_nodes, false),
          m_distances(num_nodes, std::numeric_limits<uint32_t>::max()),
          m_num_contracted_neighbors(num_nodes, 0) {}
    ContractionData() = default;
};

// define dijkstra query data for efficient memory reuse for multiple queries
struct DijkstraQueryData {
    uint32_t source;
    uint32_t target;
    uint32_t distance = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> distances;
    std::vector<uint32_t> predecessors;

    int num_pq_pops = 0;

    // keep track of the nodes that need to be reset
    // cheaper than iterating over the whole vector
    std::vector<uint32_t> reset_nodes;

    DijkstraQueryData(uint32_t num_nodes)
        : distances(num_nodes, std::numeric_limits<uint32_t>::max()), predecessors(num_nodes, -1) {}

    void reset() {
        distance = std::numeric_limits<uint32_t>::max();
        for (uint32_t node : reset_nodes) {
            distances[node] = std::numeric_limits<uint32_t>::max();
            predecessors[node] = -1;
        }
        num_pq_pops = 0;
        reset_nodes.clear();
    }

    bool needReset() { return !reset_nodes.empty(); }
};

struct BiDirectionalDijkstraQueryData {
    uint32_t source;
    uint32_t target;
    uint32_t meeting_node = -1;
    uint32_t distance = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> distances_fwd;
    std::vector<uint32_t> distances_bwd;
    std::vector<uint32_t> predecessors_fwd;
    std::vector<uint32_t> predecessors_bwd;
    std::vector<bool> visited_fwd;
    std::vector<bool> visited_bwd;

    std::vector<uint32_t> reset_nodes_fwd;
    std::vector<uint32_t> reset_nodes_bwd;

    int num_pq_pops = 0;

    BiDirectionalDijkstraQueryData(uint32_t num_nodes)
        : distances_fwd(num_nodes, std::numeric_limits<uint32_t>::max()),
          distances_bwd(num_nodes, std::numeric_limits<uint32_t>::max()),
          predecessors_fwd(num_nodes, -1),
          predecessors_bwd(num_nodes, -1),
          visited_fwd(num_nodes, false),
          visited_bwd(num_nodes, false) {}

    void reset() {
        distance = std::numeric_limits<uint32_t>::max();
        meeting_node = -1;
        for (uint32_t node : reset_nodes_fwd) {
            distances_fwd[node] = std::numeric_limits<uint32_t>::max();
            predecessors_fwd[node] = -1;
            visited_fwd[node] = false;
        }
        for (uint32_t node : reset_nodes_bwd) {
            distances_bwd[node] = std::numeric_limits<uint32_t>::max();
            predecessors_bwd[node] = -1;
            visited_bwd[node] = false;
        }
        num_pq_pops = 0;
        reset_nodes_fwd.clear();
        reset_nodes_bwd.clear();
    }
    bool needReset() { return !reset_nodes_fwd.empty() || !reset_nodes_bwd.empty(); }
};

enum Heuristic {
    IN_OUT_PRODUCT = 0,
    EDGE_DIFFERENCE = 1,
    WEIGHTED_COST = 2,
    MIXED = 3,
};
}  // namespace labosm