#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <regex>
#include <tuple>
#include <vector>

namespace labosm {

/**
 * @brief BoudingBox struct used for filtering points that are outside the water
 */
struct BoudingBox {
    double min_lat;
    double max_lat;
    double min_lon;
    double max_lon;

    BoudingBox(double min_lat, double max_lat, double min_lon, double max_lon)
        : min_lat(min_lat), max_lat(max_lat), min_lon(min_lon), max_lon(max_lon) {}

    BoudingBox() : min_lat(0), max_lat(0), min_lon(0), max_lon(0) {}

    bool contains(double lat, double lon) const {
        return lat >= min_lat && lat <= max_lat && lon >= min_lon && lon <= max_lon;
    }
};

/**
 * @brief Vec3 struct used for converting lat/lon to 3D coordinates
 */
struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double lat, double lon) {
        double lat_rad = lat * M_PI / 180.0;
        double lon_rad = lon * M_PI / 180.0;
        x = cos(lat_rad) * cos(lon_rad);
        y = cos(lat_rad) * sin(lon_rad);
        z = sin(lat_rad);
    }

    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3 normalize() const {
        double magnitude = sqrt(x * x + y * y + z * z);
        if (magnitude == 0) return Vec3(0, 0, 0);
        return Vec3(x / magnitude, y / magnitude, z / magnitude);
    }

    Vec3 negative() const { return Vec3(-x, -y, -z); }
};

constexpr double R_EARTH = 6'371'000.0;  // metres

// used for nanoflann KD-tree
struct PointXYZ {
    double x, y, z;
};

// used for nanoflann KD-tree
struct PointCloud {
    std::vector<PointXYZ> pts;

    // nanoflann API
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline double kdtree_get_pt(size_t idx, int dim) const {
        const PointXYZ& p = pts[idx];
        return dim == 0 ? p.x : (dim == 1 ? p.y : p.z);
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const {
        return false;
    }
};

// Use the harvesine formula to calculate the great-circle distance between two points on the Earth
inline int greatCircleDistance(double lat1, double lon1, double lat2, double lon2) {
    double dLat = lat2 - lat1, dLon = lon2 - lon1;
    double a = std::pow(std::sin(dLat * 0.5), 2) + std::pow(std::sin(dLon * 0.5), 2) * std::cos(lat1) * std::cos(lat2);
    double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    return static_cast<int>(R_EARTH * c);
}

// SimpleEdge struct used for the simple graph data structure for dijkstra, since we dont have underlying edges
struct SimpleEdge {
    int m_target;  // corresponds to index of node in graph data structure
    int m_cost;

    SimpleEdge(int target, int cost) : m_target(target), m_cost(cost) {}
};

struct Edge {
    int m_target;  // corresponds to index of node in graph data structure
    int m_cost;

    // represents the underlying edges if this is a shortcut
    std::tuple<int, uint16_t, bool> m_child_1;
    std::tuple<int, uint16_t, bool> m_child_2;
    Edge(int target, int cost)
        : m_target(target), m_cost(cost), m_child_1({-1, -1, false}), m_child_2({-1, -1, false}) {}
    Edge(int target, int cost, std::tuple<int, int, bool> child_1, std::tuple<int, int, bool> child_2)
        : m_target(target), m_cost(cost), m_child_1(child_1), m_child_2(child_2) {}

    bool isShortcut() const { return (std::get<0>(m_child_1) != -1) && (std::get<0>(m_child_2) != -1); }
};

// ContractionEdge struct only used when creating CH
struct ContractionEdge {
    int m_target;  // corresponds to index of node in graph data structure
    int m_cost;
    // data used to find the underlying edge if this is a shortcut
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
 * Similar to the QueryData struct, to avoid allocating/deallocating a huge amount of memory
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

/**
 * @brief Used for dijkstra distance and path queries.
 * Pre defining all data allows for many queries at a time without allocating/deallocating a huge amount of memory
 */
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

// Function to parse a line and extract a number
inline int parseLine(std::string line) {
    return stoi(std::regex_replace(line, std::regex("[^0-9]*([0-9]+).*"), std::string("$1")));
}

// Function to get current memory usage in KB
// https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
inline int getMemoryUsage() {
    int result = -1;
    std::ifstream file("/proc/self/status");
    std::string line;
    while (getline(file, line)) {
        if (line.find("VmSize") != std::string::npos) {
            result = parseLine(line);
            break;
        }
    }
    file.close();
    return result;
}
}  // namespace labosm
