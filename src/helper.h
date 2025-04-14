#pragma once

#include <vector>
#include <tuple>

namespace labosm {

struct Edge {
    int m_target;  // corresponds to index of node in graph data structure
    int m_cost;
    std::tuple<int, int, bool> m_child_1;
    std::tuple<int, int, bool> m_child_2;
    Edge(int target, int cost)
        : m_target(target), m_cost(cost), m_child_1({-1, -1, false}), m_child_2({-1, -1, false}) {}
    Edge(int target, int cost, std::tuple<int, int, bool> child_1, std::tuple<int, int, bool> child_2)
        : m_target(target), m_cost(cost), m_child_1(child_1), m_child_2(child_2) {}

    bool isShortcut() const { return (std::get<0>(m_child_1) != -1) && (std::get<0>(m_child_2) != -1); }
};

// Edge struct only used when creating CH
struct ContractionEdge {
    int m_target;  // corresponds to index of node in graph data structure
    int m_cost;
    // data used to find the underlying later on
    int m_contraction_node;
    int m_num_underlying_arcs;  // used for the mixed heuristic
    ContractionEdge(int target, int cost, int contraction_node, int num_underlying_arcs)
        : m_target(target),
          m_cost(cost),
          m_contraction_node(contraction_node),
          m_num_underlying_arcs(num_underlying_arcs) {}
    ContractionEdge(int target, int cost) : m_target(target), m_cost(cost), m_contraction_node(-1) {}

    bool isShortcut() const { return m_contraction_node != -1; }
};

/**
 * @brief Used for calculating the CH
 *
 */
struct ContractionData {
    std::vector<bool> m_outgoing;
    std::vector<int> m_reset_outgoing;
    std::vector<bool> m_visited;
    std::vector<int> m_reset_visited;
    std::vector<int> m_distances;
    std::vector<int> m_reset_distances;
    std::vector<int> m_num_contracted_neighbors;
    std::vector<std::pair<int, ContractionEdge>> m_shortcuts_fwd;
    std::vector<std::pair<int, ContractionEdge>> m_shortcuts_bwd;

    ContractionData(int num_nodes)
        : m_outgoing(num_nodes, false),
          m_visited(num_nodes, false),
          m_distances(num_nodes, std::numeric_limits<int>::max()) {}
    ContractionData() = default;
};

// define dijkstra query data for efficient memory reuse for multiple queries
struct DijkstraQueryData {
    int m_start;
    int m_end;
    int m_distance = std::numeric_limits<int>::max();
    std::vector<int> m_distances;
    std::vector<int> m_prev;
    std::vector<int> m_path;

    int num_pq_pops = 0;

    // keep track of the nodes that need to be reset
    // cheaper than iterating over the whole vector
    std::vector<int> m_reset_nodes;

    DijkstraQueryData(int num_nodes) : m_distances(num_nodes, std::numeric_limits<int>::max()), m_prev(num_nodes, -1) {}

    void reset() {
        m_distance = std::numeric_limits<int>::max();
        for (int node : m_reset_nodes) {
            m_distances[node] = std::numeric_limits<int>::max();
            m_prev[node] = -1;
        }
        num_pq_pops = 0;
        m_reset_nodes.clear();
        m_path.clear();
    }

    bool needReset() { return !m_reset_nodes.empty(); }
};

/**
 * @brief Used for distances and path queries.
 * Pre defining all data allows for many queries at a time without allocating/deallocating a huge amount of memory
 *
 */
struct QueryData {
    int m_start;
    int m_end;
    int m_meeting_node = -1;
    int m_distance = std::numeric_limits<int>::max();
    std::vector<int> m_distances_fwd;
    std::vector<int> m_distances_bwd;

    std::vector<std::pair<int, int>> m_fwd_prev_edge;
    std::vector<std::pair<int, int>> m_bwd_prev_edge;
    std::vector<int> m_shortest_path;

    std::vector<bool> m_visited_fwd;
    std::vector<bool> m_visited_bwd;

    std::vector<int> m_reset_nodes_fwd;
    std::vector<int> m_reset_nodes_bwd;

    int num_pq_pops = 0;

    QueryData(int num_nodes)
        : m_distances_fwd(num_nodes, std::numeric_limits<int>::max()),
          m_distances_bwd(num_nodes, std::numeric_limits<int>::max()),
          m_fwd_prev_edge(num_nodes),
          m_bwd_prev_edge(num_nodes),
          m_visited_fwd(num_nodes, false),
          m_visited_bwd(num_nodes, false) {}

    void reset() {
        m_distance = std::numeric_limits<int>::max();
        m_meeting_node = -1;
        for (int node : m_reset_nodes_fwd) {
            m_distances_fwd[node] = std::numeric_limits<int>::max();
            m_fwd_prev_edge[node] = {-1, -1};
            m_visited_fwd[node] = false;
        }
        for (int node : m_reset_nodes_bwd) {
            m_distances_bwd[node] = std::numeric_limits<int>::max();
            m_bwd_prev_edge[node] = {-1, -1};
            m_visited_bwd[node] = false;
        }
        num_pq_pops = 0;
        m_reset_nodes_fwd.clear();
        m_reset_nodes_bwd.clear();
        m_shortest_path.clear();
    }
    bool needReset() const { return !m_reset_nodes_fwd.empty() || !m_reset_nodes_bwd.empty(); }
};

/**
 * @brief Defines which heuristic is used for calculating contraction hierarchies
 * IN_OUT and EDGE_DIFFERENCE are baseline heuristics.
 * WEIGHTED_COST is the preferred heuristic when not using IS.
 * MIXED is preferred when using IS.
 */
enum Heuristic { IN_OUT = 0, EDGE_DIFFERENCE = 1, WEIGHTED_COST = 2, MIXED = 3 };
}  // namespace labosm
