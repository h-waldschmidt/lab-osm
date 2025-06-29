#include "graph.h"

#include <sys/types.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "../third-party/radix_heap.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace labosm {

std::string heuristicToString(Heuristic heuristic) {
    switch (heuristic) {
        case IN_OUT:
            return "IN_OUT";
        case EDGE_DIFFERENCE:
            return "EDGE_DIFFERENCE";
        case WEIGHTED_COST:
            return "WEIGHTED_COST";
        case MIXED:
            return "MIXED";
        default:
            return "UNKNOWN";
    }
}

Graph::Graph(const std::string& path, bool ch_available, int num_threads, Heuristic ch_heuristic)
    : m_ch_available(ch_available), m_num_threads(num_threads), m_is(true) {
    omp_set_num_threads(m_num_threads);

    bool is_chfmi = path.size() >= 6 && path.substr(path.size() - 6) == ".chfmi";
    if (is_chfmi && m_ch_available) {
        readGraphFromCHFMI(path);
        m_num_nodes = m_graph_contr.size();
    } else {
        if (m_ch_available)
            readGraph(path);
        else
            readSimpleGraph(path);

        if (m_ch_available)
            m_num_nodes = m_graph_edges.empty() ? 0 : static_cast<int>(m_graph_indices.size()) - 1;
        else
            m_num_nodes = static_cast<int>(m_dijkstra_indices.size()) - 1;

        if (m_ch_available) {
            createCHwithIS(ch_heuristic);
            std::string output_filename_base = path;
            size_t pos_fmi = output_filename_base.rfind(".fmi");
            if (pos_fmi != std::string::npos) {
                output_filename_base.replace(pos_fmi, 4, "");
            }
            size_t pos_chfmi = output_filename_base.rfind(".chfmi");
            if (pos_chfmi != std::string::npos) {
                output_filename_base.replace(pos_chfmi, 6, "");
            }
            std::string output_filename = output_filename_base + "_" + heuristicToString(ch_heuristic) + ".chfmi";
            writeGraphToCHFMI(output_filename);
        }
    }

    createReverseGraph();
}

void Graph::readGraph(const std::string& path) {
    std::cout << "Started reading graph file." << "\n";

    auto begin = std::chrono::high_resolution_clock::now();
    std::ifstream infile;
    try {
        infile.open(path);
        if (!infile.good()) throw std::runtime_error("File doesn't exist!");
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        exit(1);
    }

    std::string line = "#";

    // skip the metadata which begins with #
    while (line[0] == '#') getline(infile, line);

    getline(infile, line);
    int num_nodes = std::stoi(line);
    getline(infile, line);
    int num_edges = std::stoi(line);

    m_node_level.clear();

    m_node_coords.resize(num_nodes);

    for (int i = 0; i < num_nodes; ++i) {
        getline(infile, line);
        std::stringstream ss(line);
        std::string s;

        // skip unused fields
        getline(ss, s, ' ');
        getline(ss, s, ' ');

        std::pair<double, double> coords;
        getline(ss, s, ' ');
        coords.first = stod(s);
        getline(ss, s, ' ');
        coords.second = stod(s);
        m_node_coords[i] = coords;
    }

    // Initialize CH graph structures
    m_graph_indices.clear();
    m_graph_indices.resize(num_nodes + 1, 0);
    m_graph_edges.clear();

    // First pass: count edges per node to build indices
    std::vector<int> edge_counts(num_nodes, 0);
    std::streampos edge_start_pos = infile.tellg();

    for (int i = 0; i < num_edges; ++i) {
        getline(infile, line);
        std::stringstream ss(line);

        std::string s;
        getline(ss, s, ' ');
        int src = std::stoi(s);
        edge_counts[src]++;
    }

    // Build indices array (cumulative sum)
    for (int i = 0; i < num_nodes; ++i) {
        m_graph_indices[i + 1] = m_graph_indices[i] + edge_counts[i];
    }

    // Reserve space for all edges
    m_graph_edges.resize(m_graph_indices[num_nodes]);

    // Reset to beginning of edge data
    infile.seekg(edge_start_pos);

    // Second pass: actually read and store edges
    std::vector<uint32_t> current_positions = m_graph_indices;  // Copy of indices to track current position
    for (int i = 0; i < num_edges; ++i) {
        getline(infile, line);
        std::stringstream ss(line);

        std::string s;
        getline(ss, s, ' ');
        int src = std::stoi(s);
        getline(ss, s, ' ');
        int target = std::stoi(s);
        getline(ss, s, ' ');
        int cost = std::stoi(s);

        m_graph_edges[current_positions[src]++] = Edge(target, cost);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished reading graph file. Took " << elapsed.count() << " milliseconds " << "\n";
}

void Graph::readSimpleGraph(const std::string& path) {
    std::cout << "Started reading graph file for dijkstra queries." << "\n";
    auto begin = std::chrono::high_resolution_clock::now();
    std::ifstream infile;
    try {
        infile.open(path);
        if (!infile.good()) throw std::runtime_error("File doesn't exist!");
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        exit(1);
    }

    std::string line = "#";

    // skip the metadata which begins with #
    while (line[0] == '#') getline(infile, line);

    getline(infile, line);
    int num_nodes = std::stoi(line);
    getline(infile, line);
    int num_edges = std::stoi(line);

    m_node_level.clear();

    m_node_coords.resize(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
        getline(infile, line);
        std::stringstream ss(line);
        std::string s;

        // skip unused fields
        getline(ss, s, ' ');
        getline(ss, s, ' ');

        std::pair<double, double> coords;
        getline(ss, s, ' ');
        coords.first = stod(s);
        getline(ss, s, ' ');
        coords.second = stod(s);
        m_node_coords[i] = coords;
    }

    // Initialize CH graph structures (not used in simple graph mode, but initialize to avoid issues)
    m_graph_indices.clear();
    m_graph_indices.resize(num_nodes + 1, 0);
    m_graph_edges.clear();

    // Initialize dijkstra graph structures
    m_dijkstra_indices.clear();
    m_dijkstra_indices.resize(num_nodes + 1, 0);
    m_dijkstra_edges.clear();

    // First pass: count edges per node to build indices
    std::vector<int> edge_counts(num_nodes, 0);
    std::streampos edge_start_pos = infile.tellg();

    for (int i = 0; i < num_edges; ++i) {
        getline(infile, line);
        std::stringstream ss(line);

        std::string s;
        getline(ss, s, ' ');
        int src = std::stoi(s);
        edge_counts[src]++;
    }

    // Build indices array (cumulative sum)
    for (int i = 0; i < num_nodes; ++i) {
        m_dijkstra_indices[i + 1] = m_dijkstra_indices[i] + edge_counts[i];
    }

    // Reserve space for all edges
    m_dijkstra_edges.resize(m_dijkstra_indices[num_nodes]);

    // Reset to beginning of edge data
    infile.seekg(edge_start_pos);

    // Second pass: actually read and store edges
    std::vector<uint64_t> current_positions = m_dijkstra_indices;  // Copy of indices to track current position
    for (int i = 0; i < num_edges; ++i) {
        getline(infile, line);
        std::stringstream ss(line);

        std::string s;
        getline(ss, s, ' ');
        int src = std::stoi(s);
        getline(ss, s, ' ');
        int target = std::stoi(s);
        getline(ss, s, ' ');
        int cost = std::stoi(s);

        m_dijkstra_edges[current_positions[src]++] = SimpleEdge(target, cost);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished reading graph file. Took " << elapsed.count() << " milliseconds " << "\n";
}

void Graph::createReverseGraph() {
    std::cout << "Started creating reverse graph." << "\n";

    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph_edges.empty() && m_graph_contr.empty()) {
        std::cout << "Can't create reverse graph, because graph is empty" << "\n";
        return;
    }

    m_reverse_graph_indices.clear();
    m_reverse_graph_indices.resize(m_num_nodes + 1, 0);
    m_reverse_graph_edges.clear();

    if (m_ch_available)
        createReverseGraphCH();
    else
        createReverseGraphNormal();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating reverse graph. Took " << elapsed.count() << " milliseconds " << "\n";
}

void Graph::dijkstraQuery(DijkstraQueryData& data) {
    if (data.needReset()) {
        data.reset();
    }

    radix_heap::pair_radix_heap<int, int> pq;

    pq.emplace(0, data.m_start);
    data.m_prev[data.m_start] = data.m_start;

    while (!pq.empty()) {
        int dist = pq.top_key();
        int node = pq.top_value();
        pq.pop();
        data.num_pq_pops++;

        if (data.m_distances[node] <= dist) continue;

        data.m_distances[node] = dist;
        data.m_reset_nodes.push_back(node);

        if (node == data.m_end) {
            data.m_distance = dist;
            return;
        }

        for (uint64_t i = m_dijkstra_indices[node]; i < m_dijkstra_indices[node + 1]; ++i) {
            const auto& edge = m_dijkstra_edges[i];
            int new_distance = dist + edge.m_cost;
            if (new_distance < data.m_distances[edge.m_target]) {
                pq.emplace(new_distance, edge.m_target);
                data.m_prev[edge.m_target] = node;
            }
        }
    }
}

void Graph::dijkstraExtractPath(DijkstraQueryData& data) {
    if (!(data.m_distance < std::numeric_limits<int>::max() || data.m_distance > -1)) {
        std::cout << "Can't return path for invalid data!" << "\n";
        return;
    }

    data.m_path.clear();

    int cur_node = data.m_end;
    while (cur_node != data.m_start) {
        data.m_path.push_back(cur_node);
        cur_node = data.m_prev[cur_node];
    }
    data.m_path.push_back(data.m_start);
    std::reverse(data.m_path.begin(), data.m_path.end());
}

void Graph::contractionHierarchyQuery(QueryData& data) {
    if (data.m_start < 0 || data.m_end < 0) {
        std::cout << "Invalid start or end nodes!" << "\n";
        return;
    }

    if (data.needReset()) {
        data.reset();
    }

    radix_heap::pair_radix_heap<int, int> fwd_pq;
    radix_heap::pair_radix_heap<int, int> bwd_pq;

    data.m_distances_fwd[data.m_start] = 0;
    data.m_distances_bwd[data.m_end] = 0;

    // first corresponds to distance and second is node index
    fwd_pq.emplace(0, data.m_start);
    bwd_pq.emplace(0, data.m_end);

    data.m_fwd_prev_edge[data.m_start] = std::make_pair(data.m_start, -1);
    data.m_bwd_prev_edge[data.m_end] = std::make_pair(data.m_end, -1);

    data.m_visited_fwd[data.m_start] = true;
    data.m_visited_bwd[data.m_end] = true;

    data.m_reset_nodes_fwd.push_back(data.m_start);
    data.m_reset_nodes_bwd.push_back(data.m_end);

    // bidirectional dijkstra
    while (!fwd_pq.empty() || !bwd_pq.empty()) {
        while (!fwd_pq.empty()) {
            int fwd_dist = fwd_pq.top_key();
            int fwd_node_id = fwd_pq.top_value();
            fwd_pq.pop();
            data.num_pq_pops++;

            if (data.m_visited_fwd[fwd_node_id] && fwd_dist > data.m_distances_fwd[fwd_node_id]) continue;
            if (fwd_dist > data.m_distance) break;

            // forward step
            uint64_t fwd_start_idx = m_graph_indices[fwd_node_id];
            uint64_t fwd_end_idx = m_graph_indices[fwd_node_id + 1];
            for (uint64_t i = fwd_start_idx; i < fwd_end_idx; ++i) {
                Edge& e = m_graph_edges[i];

                if (!data.m_visited_fwd[e.m_target] || data.m_distances_fwd[e.m_target] > fwd_dist + e.m_cost) {
                    if (!data.m_visited_fwd[e.m_target]) data.m_reset_nodes_fwd.push_back(e.m_target);

                    data.m_distances_fwd[e.m_target] = fwd_dist + e.m_cost;
                    fwd_pq.emplace(data.m_distances_fwd[e.m_target], e.m_target);
                    data.m_visited_fwd[e.m_target] = true;

                    data.m_fwd_prev_edge[e.m_target] = std::make_pair(fwd_node_id, static_cast<int>(i - fwd_start_idx));
                }

                if (data.m_visited_bwd[e.m_target] &&
                    data.m_distances_fwd[fwd_node_id] + e.m_cost + data.m_distances_bwd[e.m_target] < data.m_distance) {
                    data.m_distance = data.m_distances_fwd[fwd_node_id] + e.m_cost + data.m_distances_bwd[e.m_target];
                    data.m_meeting_node = e.m_target;
                }
            }

            break;
        }

        while (!bwd_pq.empty()) {
            int bwd_dist = bwd_pq.top_key();
            int bwd_node_id = bwd_pq.top_value();
            bwd_pq.pop();
            data.num_pq_pops++;

            if (data.m_visited_bwd[bwd_node_id] && bwd_dist > data.m_distances_bwd[bwd_node_id]) continue;
            if (bwd_dist > data.m_distance) break;

            // backward step
            uint64_t bwd_start_idx = m_reverse_graph_indices[bwd_node_id];
            uint64_t bwd_end_idx = m_reverse_graph_indices[bwd_node_id + 1];
            for (uint64_t i = bwd_start_idx; i < bwd_end_idx; ++i) {
                Edge& e = m_reverse_graph_edges[i];
                if (!data.m_visited_bwd[e.m_target] || data.m_distances_bwd[e.m_target] > bwd_dist + e.m_cost) {
                    if (!data.m_visited_bwd[e.m_target]) data.m_reset_nodes_bwd.push_back(e.m_target);

                    data.m_distances_bwd[e.m_target] = bwd_dist + e.m_cost;
                    bwd_pq.emplace(data.m_distances_bwd[e.m_target], e.m_target);
                    data.m_visited_bwd[e.m_target] = true;

                    data.m_bwd_prev_edge[e.m_target] = std::make_pair(bwd_node_id, static_cast<int>(i - bwd_start_idx));
                }

                if (data.m_visited_fwd[e.m_target] &&
                    data.m_distances_bwd[bwd_node_id] + e.m_cost + data.m_distances_fwd[e.m_target] < data.m_distance) {
                    data.m_distance = data.m_distances_bwd[bwd_node_id] + e.m_cost + data.m_distances_fwd[e.m_target];
                    data.m_meeting_node = e.m_target;
                }
            }
            break;
        }
    }
}

void Graph::contractionHierarchyExtractPath(QueryData& data) {
    if (data.m_distance == std::numeric_limits<int>::max() || data.m_distance == -1 || data.m_meeting_node == -1) {
        std::cout << "Can't return path for invalid data!" << "\n";
        return;
    }

    data.m_shortest_path.clear();

    // Use temporary vectors from QueryData to avoid repeated allocations
    data.m_temp_fwd_edges.clear();
    data.m_temp_bwd_edges.clear();
    data.m_temp_fwd_path.clear();
    data.m_temp_bwd_path.clear();

    int cur_node = data.m_meeting_node;

    // Gather forward edges - use temporary vector
    while (cur_node != data.m_start) {
        data.m_temp_fwd_edges.emplace_back(data.m_fwd_prev_edge[cur_node].first, data.m_fwd_prev_edge[cur_node].second,
                                           true);
        cur_node = data.m_fwd_prev_edge[cur_node].first;
    }
    std::reverse(data.m_temp_fwd_edges.begin(), data.m_temp_fwd_edges.end());

    // Unpack forward edges efficiently - reserve space to avoid reallocations
    if (data.m_temp_fwd_edges.empty()) {
        data.m_temp_fwd_path.push_back(data.m_start);
    } else {
        // Estimate path length to reserve space
        data.m_temp_fwd_path.reserve(data.m_temp_fwd_edges.size() * 3);

        unpackEdge(data.m_temp_fwd_edges[0], data.m_temp_fwd_path);

        // Use a small temporary vector for each edge unpacking to avoid frequent reallocations
        std::vector<int> temp_edge_path;
        temp_edge_path.reserve(50);  // Reserve reasonable space for typical shortcut depth

        for (size_t i = 1; i < data.m_temp_fwd_edges.size(); ++i) {
            temp_edge_path.clear();
            unpackEdge(data.m_temp_fwd_edges[i], temp_edge_path);
            // More efficient than insert: manually append elements
            data.m_temp_fwd_path.reserve(data.m_temp_fwd_path.size() + temp_edge_path.size() - 1);
            for (size_t j = 1; j < temp_edge_path.size(); ++j) {
                data.m_temp_fwd_path.push_back(temp_edge_path[j]);
            }
        }
    }

    // Gather backward edges - use temporary vector
    cur_node = data.m_meeting_node;
    while (cur_node != data.m_end) {
        data.m_temp_bwd_edges.emplace_back(data.m_bwd_prev_edge[cur_node].first, data.m_bwd_prev_edge[cur_node].second,
                                           false);
        cur_node = data.m_bwd_prev_edge[cur_node].first;
    }

    // Unpack backward edges efficiently
    if (data.m_temp_bwd_edges.empty()) {
        data.m_temp_bwd_path.push_back(data.m_end);
    } else {
        // Estimate path length to reserve space
        data.m_temp_bwd_path.reserve(data.m_temp_bwd_edges.size() * 3);

        unpackEdge(data.m_temp_bwd_edges[0], data.m_temp_bwd_path);

        // Reuse the same temporary vector
        std::vector<int> temp_edge_path;
        temp_edge_path.reserve(50);

        for (size_t i = 1; i < data.m_temp_bwd_edges.size(); ++i) {
            temp_edge_path.clear();
            unpackEdge(data.m_temp_bwd_edges[i], temp_edge_path);
            // More efficient than insert: manually append elements
            data.m_temp_bwd_path.reserve(data.m_temp_bwd_path.size() + temp_edge_path.size() - 1);
            for (size_t j = 1; j < temp_edge_path.size(); ++j) {
                data.m_temp_bwd_path.push_back(temp_edge_path[j]);
            }
        }
    }

    // Combine paths efficiently - reserve total space upfront
    data.m_shortest_path.reserve(data.m_temp_fwd_path.size() + data.m_temp_bwd_path.size() - 1);

    // Copy forward path
    for (int node : data.m_temp_fwd_path) {
        data.m_shortest_path.push_back(node);
    }

    // Copy backward path (skip first element to avoid duplication)
    for (size_t i = 1; i < data.m_temp_bwd_path.size(); ++i) {
        data.m_shortest_path.push_back(data.m_temp_bwd_path[i]);
    }
}

void Graph::unpackEdge(const std::tuple<int, int, bool>& edge_index, std::vector<int>& path) {
    const Edge* edge = nullptr;
    // the bool indicates if we need to look at bwd or fwd edge/graph
    if (std::get<2>(edge_index)) {
        // Forward graph: calculate absolute index from node and relative edge index
        uint64_t node_start = m_graph_indices[std::get<0>(edge_index)];
        edge = &m_graph_edges[node_start + std::get<1>(edge_index)];
    } else {
        // Reverse graph: calculate absolute index from node and relative edge index
        uint64_t node_start = m_reverse_graph_indices[std::get<0>(edge_index)];
        edge = &m_reverse_graph_edges[node_start + std::get<1>(edge_index)];
    }

    if (!edge->isShortcut()) {
        if (path.empty() || (path.back() != std::get<0>(edge_index) && path.back() != edge->m_target)) {
            if (std::get<2>(edge_index))
                path.push_back(std::get<0>(edge_index));
            else
                path.push_back(edge->m_target);
        }

        if (std::get<2>(edge_index))
            path.push_back(edge->m_target);
        else
            path.push_back(std::get<0>(edge_index));

    } else {
        if (std::get<2>(edge_index)) {
            unpackEdge(edge->m_child_1, path);
            unpackEdge(edge->m_child_2, path);
        } else {
            unpackEdge(edge->m_child_2, path);
            unpackEdge(edge->m_child_1, path);
        }
    }
}

void Graph::createHubLabelsWithIS() {
    std::cout << "Started creating hub labels." << "\n";
    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph_edges.empty() || m_reverse_graph_edges.empty()) {
        std::cout << "Can't create hub labels, because graph or reverse graph is empty" << "\n";
        return;
    }

    m_fwd_indices.clear();
    m_fwd_indices.resize(m_num_nodes + 1, 0);
    m_bwd_indices.clear();
    m_bwd_indices.resize(m_num_nodes + 1, 0);

    m_fwd_hub_labels.clear();
    m_bwd_hub_labels.clear();

    // Reserve space to avoid frequent reallocations during construction
    // Estimate based on typical hub label density
    const size_t estimated_total_labels = static_cast<size_t>(m_num_nodes) * 300;  // Rough estimate
    m_fwd_hub_labels.reserve(estimated_total_labels);
    m_bwd_hub_labels.reserve(estimated_total_labels);

    // Since nodes are already reordered by importance in createReverseGraphCH
    // (node 0 = highest importance), we can process them in order without additional sorting

    // Group nodes by level to process level by level (for parallel efficiency)
    std::unordered_map<int, std::vector<int>> level_buckets;
    for (int i = 0; i < m_num_nodes; i++) {
        level_buckets[m_node_level[i]].push_back(i);
    }

    // Extract sorted levels (highest first)
    std::vector<int> different_levels_vec;
    different_levels_vec.reserve(level_buckets.size());
    for (const auto& pair : level_buckets) {
        different_levels_vec.push_back(pair.first);
    }
    std::sort(different_levels_vec.begin(), different_levels_vec.end(), std::greater<int>());

    // Convert to vector of vectors for processing
    std::vector<std::vector<int>> level_buckets_vec;
    level_buckets_vec.reserve(different_levels_vec.size());
    for (int level : different_levels_vec) {
        level_buckets_vec.push_back(std::move(level_buckets[level]));
    }

    uint32_t num_calculated = 0;

    const int MAX_THREADS = 16;

    // reserve space for fwd_labels/bwd_labels to avoid alocations for each node
    std::vector<std::vector<std::tuple<uint32_t, uint32_t, uint16_t, uint16_t>>> fwd_labels_vec(MAX_THREADS);
    std::vector<std::vector<std::tuple<uint32_t, uint32_t, uint16_t, uint16_t>>> bwd_labels_vec(MAX_THREADS);

    for (int i = 0; i < MAX_THREADS; ++i) {
        fwd_labels_vec[i].reserve(300);
        bwd_labels_vec[i].reserve(300);
    }

    // at the beginning set threads to 1
    omp_set_num_threads(1);

    for (int i = 0; i < different_levels_vec.size(); ++i) {
        if (i == 0) {
            // Handle the highest level node (should be node 0 due to reordering)
            int node = level_buckets_vec[i][0];
            m_fwd_hub_labels.emplace_back(static_cast<uint32_t>(node), 0U, static_cast<uint16_t>(-1),
                                          static_cast<uint32_t>(-1));
            m_bwd_hub_labels.emplace_back(static_cast<uint32_t>(node), 0U, static_cast<uint16_t>(-1),
                                          static_cast<uint32_t>(-1));
            // Since nodes are ordered by importance, node indices equal node IDs
            m_fwd_indices[node + 1] = 1;
            m_bwd_indices[node + 1] = 1;
            continue;
        }

#pragma omp parallel for ordered schedule(static, 1)
        for (int j = 0; j < static_cast<int>(level_buckets_vec[i].size()); ++j) {
            int node = level_buckets_vec[i][j];

            auto& fwd_labels = fwd_labels_vec[omp_get_thread_num()];
            auto& bwd_labels = bwd_labels_vec[omp_get_thread_num()];
            fwd_labels.clear();
            bwd_labels.clear();

            fwd_labels.emplace_back(static_cast<uint32_t>(node), 0U, static_cast<uint16_t>(-1),
                                    static_cast<uint32_t>(-1));
            bwd_labels.emplace_back(static_cast<uint32_t>(node), 0U, static_cast<uint16_t>(-1),
                                    static_cast<uint32_t>(-1));

            // fwd labels - cache optimization to reduce repeated lookups
            const int node_level = m_node_level[node];
            const uint64_t graph_start_idx = m_graph_indices[node];
            const uint64_t graph_end_idx = m_graph_indices[node + 1];

            for (uint64_t k = graph_start_idx; k < graph_end_idx; ++k) {
                const Edge& edge = m_graph_edges[k];

                // Since nodes are ordered by importance, target node index equals target node ID
                const int target_node = edge.m_target;
                for (uint64_t j = m_fwd_indices[target_node]; j < m_fwd_indices[target_node + 1]; j++) {
                    fwd_labels.emplace_back(std::get<0>(m_fwd_hub_labels[j]),
                                            std::get<1>(m_fwd_hub_labels[j]) + static_cast<uint32_t>(edge.m_cost),
                                            static_cast<uint16_t>(k - graph_start_idx),  // Store relative edge index
                                            static_cast<uint16_t>(j - m_fwd_indices[target_node]));
                }
            }

            // remove duplicates - optimized version
            std::sort(fwd_labels.begin(), fwd_labels.end(), [](const auto& left, const auto& right) {
                if (std::get<0>(left) != std::get<0>(right)) return std::get<0>(left) < std::get<0>(right);
                return std::get<1>(left) < std::get<1>(right);
            });

            // Remove duplicates by keeping only the best (smallest distance) for each node
            auto write_it = fwd_labels.begin();
            for (auto read_it = fwd_labels.begin(); read_it != fwd_labels.end();) {
                *write_it = std::move(*read_it);
                auto next_it = read_it + 1;
                // Skip all duplicates with worse distances
                while (next_it != fwd_labels.end() && std::get<0>(*next_it) == std::get<0>(*write_it)) {
                    if (std::get<1>(*next_it) < std::get<1>(*write_it)) {
                        *write_it = std::move(*next_it);
                    }
                    ++next_it;
                }
                read_it = next_it;
                ++write_it;
            }
            fwd_labels.erase(write_it, fwd_labels.end());

            for (auto iter = fwd_labels.begin(); iter != fwd_labels.end();) {
                int best_dist = std::numeric_limits<int>::max();
                if (static_cast<uint32_t>(node) != std::get<0>(*iter))
                    best_dist = simplifiedHubLabelQuery(fwd_labels, static_cast<int>(std::get<0>(*iter)));
                if (best_dist < std::get<1>(*iter))
                    iter = fwd_labels.erase(iter);
                else
                    ++iter;
            }

            // bwd labels - cache optimization to reduce repeated lookups
            const uint64_t reverse_graph_start_idx = m_reverse_graph_indices[node];
            const uint64_t reverse_graph_end_idx = m_reverse_graph_indices[node + 1];

            for (uint64_t k = reverse_graph_start_idx; k < reverse_graph_end_idx; ++k) {
                const Edge& edge = m_reverse_graph_edges[k];

                // Since nodes are ordered by importance, target node index equals target node ID
                const int target_node = edge.m_target;
                for (uint64_t j = m_bwd_indices[target_node]; j < m_bwd_indices[target_node + 1]; j++) {
                    bwd_labels.emplace_back(
                        std::get<0>(m_bwd_hub_labels[j]),
                        std::get<1>(m_bwd_hub_labels[j]) + static_cast<uint32_t>(edge.m_cost),
                        static_cast<uint16_t>(k - reverse_graph_start_idx),  // Store relative edge index
                        static_cast<uint16_t>(j - m_bwd_indices[target_node]));
                }
            }

            // remove duplicates - optimized version
            std::sort(bwd_labels.begin(), bwd_labels.end(), [](const auto& left, const auto& right) {
                if (std::get<0>(left) != std::get<0>(right)) return std::get<0>(left) < std::get<0>(right);
                return std::get<1>(left) < std::get<1>(right);
            });

            // Remove duplicates by keeping only the best (smallest distance) for each node
            write_it = bwd_labels.begin();
            for (auto read_it = bwd_labels.begin(); read_it != bwd_labels.end();) {
                *write_it = std::move(*read_it);
                auto next_it = read_it + 1;
                // Skip all duplicates with worse distances
                while (next_it != bwd_labels.end() && std::get<0>(*next_it) == std::get<0>(*write_it)) {
                    if (std::get<1>(*next_it) < std::get<1>(*write_it)) {
                        *write_it = std::move(*next_it);
                    }
                    ++next_it;
                }
                read_it = next_it;
                ++write_it;
            }
            bwd_labels.erase(write_it, bwd_labels.end());

            for (auto iter = bwd_labels.begin(); iter != bwd_labels.end();) {
                int best_dist = std::numeric_limits<int>::max();
                if (static_cast<uint32_t>(node) != std::get<0>(*iter))
                    best_dist = simplifiedHubLabelQuery(static_cast<int>(std::get<0>(*iter)), bwd_labels);
                if (best_dist < std::get<1>(*iter))
                    iter = bwd_labels.erase(iter);
                else
                    ++iter;
            }

#pragma omp ordered
            {
                // Since nodes are already ordered by importance, node index equals node ID
                m_fwd_indices[node + 1] = m_fwd_indices[node] + static_cast<uint64_t>(fwd_labels.size());
                m_bwd_indices[node + 1] = m_bwd_indices[node] + static_cast<uint64_t>(bwd_labels.size());

                for (const auto& label : fwd_labels) {
                    m_fwd_hub_labels.emplace_back(label);
                }
                for (const auto& label : bwd_labels) {
                    m_bwd_hub_labels.emplace_back(label);
                }
            }
        }
        num_calculated += static_cast<uint32_t>(level_buckets_vec[i].size());
        std::cout << "Finished Hub Labels for " << num_calculated << " num of nodes." << "\n";

        const double PARALLEL_THRESHOLD = 0.05;
        if (num_calculated > PARALLEL_THRESHOLD * m_num_nodes) {
            omp_set_num_threads(MAX_THREADS);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating hub labels. Took " << elapsed.count() << " milliseconds " << "\n";
}

std::pair<uint64_t, uint64_t> Graph::hubLabelQuery(QueryData& data) {
    if (data.m_start < 0 || data.m_end < 0) {
        std::cout << "Invalid start or end nodes!" << "\n";
        return std::make_pair(-1, -1);
    }

    if (data.needReset()) {
        data.reset();
    }

    std::cout << "Started Hub Query." << "\n";

    auto begin = std::chrono::high_resolution_clock::now();

    data.m_distance = std::numeric_limits<int>::max();
    data.m_meeting_node = -1;

    uint64_t fwd_node_index = m_fwd_indices[data.m_start];
    uint64_t fwd_next_index = m_fwd_indices[data.m_start + 1];

    uint64_t bwd_node_index = m_bwd_indices[data.m_end];
    uint64_t bwd_next_index = m_bwd_indices[data.m_end + 1];

    uint64_t best_fwd_index = fwd_node_index;
    uint64_t best_bwd_index = bwd_node_index;
    while (fwd_node_index < fwd_next_index && bwd_node_index < bwd_next_index) {
        if (std::get<0>(m_fwd_hub_labels[fwd_node_index]) == std::get<0>(m_bwd_hub_labels[bwd_node_index])) {
            if (std::get<1>(m_fwd_hub_labels[fwd_node_index]) + std::get<1>(m_bwd_hub_labels[bwd_node_index]) <
                data.m_distance) {
                data.m_meeting_node = static_cast<int>(std::get<0>(m_fwd_hub_labels[fwd_node_index]));
                data.m_distance = static_cast<int>(std::get<1>(m_fwd_hub_labels[fwd_node_index]) +
                                                   std::get<1>(m_bwd_hub_labels[bwd_node_index]));

                best_fwd_index = fwd_node_index;
                best_bwd_index = bwd_node_index;
            }
            ++fwd_node_index;
            ++bwd_node_index;
        } else if (std::get<0>(m_fwd_hub_labels[fwd_node_index]) < std::get<0>(m_bwd_hub_labels[bwd_node_index])) {
            ++fwd_node_index;
        } else {
            ++bwd_node_index;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "Finished Hub Label query. Took " << elapsed.count() << " nanoseconds " << "\n";

    return std::make_pair(best_fwd_index, best_bwd_index);
}

// NOTE: Is very similar to CH Path extraction
// But need to look at the fwd and bwd labels
void Graph::hubLabelExtractPath(QueryData& data, std::pair<uint64_t, uint64_t> hub_indices) {
    if (data.m_distance == std::numeric_limits<int>::max() || data.m_distance == -1 || data.m_meeting_node == -1) {
        std::cout << "Can't return path for invalid data!" << "\n";
        return;
    }

    data.m_shortest_path.clear();

    // Use temporary vectors from QueryData to avoid repeated allocations
    data.m_temp_fwd_edges.clear();
    data.m_temp_bwd_edges.clear();
    data.m_temp_fwd_path.clear();
    data.m_temp_bwd_path.clear();

    // Gather forward edges - use temporary vector
    int cur_node = data.m_start;
    auto cur_fwd_label = m_fwd_hub_labels[hub_indices.first];
    while (cur_node != data.m_meeting_node) {
        // Get edge using flat structure: base index + relative edge index
        uint64_t absolute_edge_idx = m_graph_indices[cur_node] + std::get<2>(cur_fwd_label);
        Edge e = m_graph_edges[absolute_edge_idx];
        data.m_temp_fwd_edges.emplace_back(cur_node, std::get<2>(cur_fwd_label), true);
        cur_node = e.m_target;
        // Calculate absolute index: base index for target node + offset stored in tuple
        cur_fwd_label = m_fwd_hub_labels[m_fwd_indices[cur_node] + static_cast<uint64_t>(std::get<3>(cur_fwd_label))];
    }

    // Unpack forward edges efficiently - reserve space to avoid reallocations
    if (data.m_temp_fwd_edges.empty()) {
        data.m_temp_fwd_path.push_back(data.m_start);
    } else {
        // Estimate path length to reserve space
        data.m_temp_fwd_path.reserve(data.m_temp_fwd_edges.size() * 3);

        unpackEdge(data.m_temp_fwd_edges[0], data.m_temp_fwd_path);

        // Use a small temporary vector for each edge unpacking to avoid frequent reallocations
        std::vector<int> temp_edge_path;
        temp_edge_path.reserve(50);  // Reserve reasonable space for typical shortcut depth

        for (size_t i = 1; i < data.m_temp_fwd_edges.size(); ++i) {
            temp_edge_path.clear();
            unpackEdge(data.m_temp_fwd_edges[i], temp_edge_path);
            // More efficient than insert: manually append elements
            data.m_temp_fwd_path.reserve(data.m_temp_fwd_path.size() + temp_edge_path.size() - 1);
            for (size_t j = 1; j < temp_edge_path.size(); ++j) {
                data.m_temp_fwd_path.push_back(temp_edge_path[j]);
            }
        }
    }

    // Gather backward edges - use temporary vector
    cur_node = data.m_end;
    auto cur_bwd_label = m_bwd_hub_labels[hub_indices.second];
    while (cur_node != data.m_meeting_node) {
        // Get edge using flat structure: base index + relative edge index
        uint64_t absolute_edge_idx = m_reverse_graph_indices[cur_node] + std::get<2>(cur_bwd_label);
        Edge e = m_reverse_graph_edges[absolute_edge_idx];
        data.m_temp_bwd_edges.emplace_back(cur_node, std::get<2>(cur_bwd_label), false);
        cur_node = e.m_target;
        // Calculate absolute index: base index for target node + offset stored in tuple
        cur_bwd_label = m_bwd_hub_labels[m_bwd_indices[cur_node] + static_cast<uint64_t>(std::get<3>(cur_bwd_label))];
    }

    std::reverse(data.m_temp_bwd_edges.begin(), data.m_temp_bwd_edges.end());

    // Unpack backward edges efficiently
    if (data.m_temp_bwd_edges.empty()) {
        data.m_temp_bwd_path.push_back(data.m_end);
    } else {
        // Estimate path length to reserve space
        data.m_temp_bwd_path.reserve(data.m_temp_bwd_edges.size() * 3);

        unpackEdge(data.m_temp_bwd_edges[0], data.m_temp_bwd_path);

        // Reuse the same temporary vector
        std::vector<int> temp_edge_path;
        temp_edge_path.reserve(50);

        for (size_t i = 1; i < data.m_temp_bwd_edges.size(); ++i) {
            temp_edge_path.clear();
            unpackEdge(data.m_temp_bwd_edges[i], temp_edge_path);
            // More efficient than insert: manually append elements
            data.m_temp_bwd_path.reserve(data.m_temp_bwd_path.size() + temp_edge_path.size() - 1);
            for (size_t j = 1; j < temp_edge_path.size(); ++j) {
                data.m_temp_bwd_path.push_back(temp_edge_path[j]);
            }
        }
    }

    // Combine paths efficiently - reserve total space upfront
    data.m_shortest_path.reserve(data.m_temp_fwd_path.size() + data.m_temp_bwd_path.size() - 1);

    // Copy forward path
    for (int node : data.m_temp_fwd_path) {
        data.m_shortest_path.push_back(node);
    }

    // Copy backward path (skip first element to avoid duplication)
    for (size_t i = 1; i < data.m_temp_bwd_path.size(); ++i) {
        data.m_shortest_path.push_back(data.m_temp_bwd_path[i]);
    }
}

int Graph::greatCircleDistance(double lat_1, double lon_1, double lat_2, double lon_2) {
    // degrees to radians
    lat_1 *= (M_PI / 180.0);
    lon_1 *= (M_PI / 180.0);
    lat_2 *= (M_PI / 180.0);
    lon_2 *= (M_PI / 180.0);

    // Using Haversine Distance: https://en.wikipedia.org/wiki/Haversine_formula
    // example in JS: https://github.com/njj/haversine/blob/develop/haversine.js
    double d_lat = lat_2 - lat_1;
    double d_lon = lon_2 - lon_1;
    double a = pow(sin(d_lat / 2.0), 2) + pow(sin(d_lon / 2.0), 2) * cos(lat_1) * cos(lat_2);
    double c = 2.0 * atan2(sqrt(a), sqrt(1 - a));

    return 6371000 * c;  // return in meters
}

double Graph::averageLabelSize() {
    double avg_label_size = 0.0;

    for (int i = 0; i < m_num_nodes; i++) {
        avg_label_size += m_fwd_indices[i + 1] - m_fwd_indices[i];
        avg_label_size += m_bwd_indices[i + 1] - m_bwd_indices[i];
    }
    return avg_label_size / static_cast<double>(2 * m_num_nodes);
}

int Graph::maxLabelSize() {
    int max_label_size = 0;

    for (int i = 0; i < m_num_nodes; i++) {
        if (m_fwd_indices[i + 1] - m_fwd_indices[i] > max_label_size)
            max_label_size = m_fwd_indices[i + 1] - m_fwd_indices[i];
        if (m_bwd_indices[i + 1] - m_bwd_indices[i] > max_label_size)
            max_label_size = m_bwd_indices[i + 1] - m_bwd_indices[i];
    }
    return max_label_size;
}

int Graph::getNearestNode(double latitude, double longitude) {
    int nearest_node = -1;
    double min_distance = std::numeric_limits<double>::max();

    // just brute force it
    // TODO: Although could reuse the kd-tree from the graph_creator
    for (int i = 0; i < m_num_nodes; ++i) {
        double distance =
            Graph::greatCircleDistance(latitude, longitude, m_node_coords[i].first, m_node_coords[i].second);
        if (distance < min_distance) {
            min_distance = distance;
            nearest_node = i;
        }
    }
    return nearest_node;
}

void Graph::createReverseGraphNormal() {
    // First pass: count edges per node for reverse graph
    std::vector<int> reverse_edge_counts(m_num_nodes, 0);

    for (int i = 0; i < m_num_nodes; ++i) {
        uint64_t start_idx = m_graph_indices[i];
        uint64_t end_idx = m_graph_indices[i + 1];
        for (uint64_t j = start_idx; j < end_idx; ++j) {
            const Edge& e = m_graph_edges[j];
            reverse_edge_counts[e.m_target]++;
        }
    }

    // Build reverse graph indices
    for (int i = 0; i < m_num_nodes; ++i) {
        m_reverse_graph_indices[i + 1] = m_reverse_graph_indices[i] + reverse_edge_counts[i];
    }

    // Reserve space for reverse edges
    m_reverse_graph_edges.resize(m_reverse_graph_indices[m_num_nodes]);

    // Second pass: populate reverse graph
    std::vector<uint32_t> current_positions = m_reverse_graph_indices;
    for (int i = 0; i < m_num_nodes; ++i) {
        uint32_t start_idx = m_graph_indices[i];
        uint32_t end_idx = m_graph_indices[i + 1];
        for (uint32_t j = start_idx; j < end_idx; ++j) {
            const Edge& e = m_graph_edges[j];
            int new_source = e.m_target;
            Edge reverse_edge = e;
            reverse_edge.m_target = i;

            m_reverse_graph_edges[current_positions[new_source]++] = reverse_edge;
        }
    }
}

void Graph::createReverseGraphCH() {
    // Create mapping from old node IDs to new node IDs based on importance (highest importance gets ID 0)
    std::vector<std::pair<int, int>> node_level_pairs;
    for (int i = 0; i < m_num_nodes; ++i) {
        node_level_pairs.push_back({m_node_level[i], i});
    }

    // Sort by level (highest level first), then by original node ID for stable ordering
    std::sort(node_level_pairs.begin(), node_level_pairs.end(),
              [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                  if (a.first != b.first) return a.first > b.first;  // Higher level first
                  return a.second < b.second;                        // Stable ordering for same level
              });

    // Create mapping arrays
    m_old_to_new_mapping.resize(m_num_nodes);          // m_old_to_new_mapping[old_id] = new_id
    std::vector<int> new_to_old_mapping(m_num_nodes);  // new_to_old_mapping[new_id] = old_id (local variable)

    for (int i = 0; i < m_num_nodes; ++i) {
        int old_id = node_level_pairs[i].second;
        int new_id = i;
        m_old_to_new_mapping[old_id] = new_id;
        new_to_old_mapping[new_id] = old_id;
    }

    // Create reordered copies of the contraction graphs with new node IDs
    std::vector<std::vector<ContractionEdge>> reordered_graph_contr(m_num_nodes);

    for (int old_i = 0; old_i < m_num_nodes; ++old_i) {
        int new_i = m_old_to_new_mapping[old_i];

        // Copy all edges with targets and contraction nodes remapped to new IDs
        for (const auto& edge : m_graph_contr[old_i]) {
            ContractionEdge new_edge = edge;
            new_edge.m_target = m_old_to_new_mapping[edge.m_target];
            // Update contraction node to new ID if it's a shortcut
            if (new_edge.isShortcut()) {
                new_edge.m_contraction_node = m_old_to_new_mapping[new_edge.m_contraction_node];
            }
            reordered_graph_contr[new_i].push_back(new_edge);
        }
    }

    // Replace the original structure with reordered one
    m_graph_contr = std::move(reordered_graph_contr);

    // Update node_level array to match new ordering BEFORE edge separation
    std::vector<int> new_node_level(m_num_nodes);
    for (int new_i = 0; new_i < m_num_nodes; ++new_i) {
        int old_i = new_to_old_mapping[new_i];
        new_node_level[new_i] = m_node_level[old_i];
    }
    m_node_level = std::move(new_node_level);

    // Now separate upward/downward edges using the NEW node ordering
    m_reverse_graph_contr.clear();
    m_reverse_graph_contr.resize(m_num_nodes);

    for (int new_i = 0; new_i < m_num_nodes; ++new_i) {
        for (auto iter = m_graph_contr[new_i].begin(); iter != m_graph_contr[new_i].end();) {
            if (m_node_level[new_i] > m_node_level[iter->m_target]) {
                ContractionEdge edge = *iter;
                int new_source = edge.m_target;
                edge.m_target = new_i;
                m_reverse_graph_contr[new_source].push_back(edge);
                iter = m_graph_contr[new_i].erase(iter);
            } else {
                iter++;
            }
        }
    }

    // Now reorder everything using the new node mapping
    m_graph_indices.clear();
    m_graph_indices.resize(m_num_nodes + 1, 0);
    m_graph_edges.clear();
    m_reverse_graph_indices.clear();
    m_reverse_graph_indices.resize(m_num_nodes + 1, 0);

    // Update node coordinates to match new ordering
    if (!m_node_coords.empty()) {
        std::vector<std::pair<double, double>> new_node_coords(m_num_nodes);
        for (int new_i = 0; new_i < m_num_nodes; ++new_i) {
            int old_i = new_to_old_mapping[new_i];
            new_node_coords[new_i] = m_node_coords[old_i];
        }
        m_node_coords = std::move(new_node_coords);
    }

    // First pass: count edges per node
    std::vector<int> fwd_edge_counts(m_num_nodes, 0);
    std::vector<int> bwd_edge_counts(m_num_nodes, 0);

    for (int i = 0; i < m_num_nodes; ++i) {
        fwd_edge_counts[i] = m_graph_contr[i].size();
        bwd_edge_counts[i] = m_reverse_graph_contr[i].size();
    }

    // Build indices arrays
    for (int i = 0; i < m_num_nodes; ++i) {
        m_graph_indices[i + 1] = m_graph_indices[i] + fwd_edge_counts[i];
        m_reverse_graph_indices[i + 1] = m_reverse_graph_indices[i] + bwd_edge_counts[i];
    }

    // Reserve space for all edges
    m_graph_edges.resize(m_graph_indices[m_num_nodes]);
    m_reverse_graph_edges.resize(m_reverse_graph_indices[m_num_nodes]);

    // Copy the data to flat structures
    for (int i = 0; i < m_num_nodes; ++i) {
        uint64_t fwd_start = m_graph_indices[i];
        uint64_t bwd_start = m_reverse_graph_indices[i];

        for (size_t j = 0; j < m_graph_contr[i].size(); ++j) {
            const auto& edge = m_graph_contr[i][j];
            m_graph_edges[fwd_start + j] = Edge(edge.m_target, edge.m_cost);
        }

        for (size_t j = 0; j < m_reverse_graph_contr[i].size(); ++j) {
            const auto& edge = m_reverse_graph_contr[i][j];
            m_reverse_graph_edges[bwd_start + j] = Edge(edge.m_target, edge.m_cost);
        }
    }

    // find the correct shortcut pointers
    for (int i = 0; i < m_num_nodes; ++i) {
        uint64_t fwd_start = m_graph_indices[i];

        for (size_t j = 0; j < m_graph_contr[i].size(); ++j) {
            const ContractionEdge& edge = m_graph_contr[i][j];
            if (edge.isShortcut()) {
                // child 1 must be in the downward graph, while child 2 must be in the upward graph

                std::tuple<int, int, bool> child_1;
                uint64_t rev_start = m_reverse_graph_indices[edge.m_contraction_node];
                uint64_t rev_end = m_reverse_graph_indices[edge.m_contraction_node + 1];
                for (uint64_t k = rev_start; k < rev_end; ++k) {
                    if (m_reverse_graph_edges[k].m_target == i) {
                        child_1 = std::make_tuple(edge.m_contraction_node, static_cast<uint16_t>(k - rev_start), false);
                        break;
                    }
                }

                std::tuple<int, int, bool> child_2;
                uint64_t fwd_start_child = m_graph_indices[edge.m_contraction_node];
                uint64_t fwd_end_child = m_graph_indices[edge.m_contraction_node + 1];
                for (uint64_t k = fwd_start_child; k < fwd_end_child; ++k) {
                    if (m_graph_edges[k].m_target == edge.m_target) {
                        child_2 =
                            std::make_tuple(edge.m_contraction_node, static_cast<uint16_t>(k - fwd_start_child), true);
                        break;
                    }
                }
                m_graph_edges[fwd_start + j].m_child_1 = child_1;
                m_graph_edges[fwd_start + j].m_child_2 = child_2;
            }
        }

        uint64_t bwd_start = m_reverse_graph_indices[i];
        for (size_t j = 0; j < m_reverse_graph_contr[i].size(); ++j) {
            const ContractionEdge& edge = m_reverse_graph_contr[i][j];
            if (edge.isShortcut()) {
                // child 1 must be in the upward graph, while child 2 must be in the downward graph
                std::tuple<int, int, bool> child_1;
                uint64_t fwd_start_child = m_graph_indices[edge.m_contraction_node];
                uint64_t fwd_end_child = m_graph_indices[edge.m_contraction_node + 1];
                for (uint64_t k = fwd_start_child; k < fwd_end_child; ++k) {
                    if (m_graph_edges[k].m_target == i) {
                        child_1 =
                            std::make_tuple(edge.m_contraction_node, static_cast<uint16_t>(k - fwd_start_child), true);
                        break;
                    }
                }

                std::tuple<int, int, bool> child_2;
                uint64_t rev_start_child = m_reverse_graph_indices[edge.m_contraction_node];
                uint64_t rev_end_child = m_reverse_graph_indices[edge.m_contraction_node + 1];
                for (uint64_t k = rev_start_child; k < rev_end_child; ++k) {
                    if (m_reverse_graph_edges[k].m_target == edge.m_target) {
                        child_2 =
                            std::make_tuple(edge.m_contraction_node, static_cast<uint16_t>(k - rev_start_child), false);
                        break;
                    }
                }
                m_reverse_graph_edges[bwd_start + j].m_child_1 = child_1;
                m_reverse_graph_edges[bwd_start + j].m_child_2 = child_2;
            }
        }
    }

    m_graph_contr.clear();
    m_reverse_graph_contr.clear();
}

void Graph::createCHwithIS(Heuristic heuristic) {
    std::cout << "Started creating CH." << "\n";

    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph_edges.empty()) {
        std::cout << "Can't create CH, because graph is empty." << "\n";
        return;
    }
    m_node_level.clear();
    m_node_level.resize(m_num_nodes);

    createContractionGraphs();

    std::vector<bool> contracted(m_num_nodes, false);
    int num_contracted = 0;
    int cur_level = 0;

    m_contr_data.resize(m_num_threads);
    for (int i = 0; i < m_num_threads; ++i) {
        m_contr_data[i] = ContractionData(m_num_nodes);
        m_contr_data[i].m_num_contracted_neighbors.resize(m_num_nodes);
    }

    // initialize importance for all nodes
    radix_heap::pair_radix_heap<int, int> importance_pq;
    for (int i = 0; i < m_num_nodes; ++i) {
        int importance = 0;
        switch (heuristic) {
            case Heuristic::IN_OUT:
                importance = inOutProductHeuristic(contracted, i);
                break;
            case Heuristic::EDGE_DIFFERENCE:
                importance = edgeDifferenceHeuristic(contracted, i);
                break;
            case Heuristic::WEIGHTED_COST:
                importance = weightedCostHeuristic(contracted, i);
                break;
            case Heuristic::MIXED:
                importance = mixedHeuristic(contracted, i, 0);
                break;
            default:
                break;
        }
        importance_pq.emplace(importance, i);
    }

    while (num_contracted != m_num_nodes) {
        // create independet set
        std::vector<bool> marked(m_num_nodes, false);
        std::vector<int> independent_set;
        radix_heap::pair_radix_heap<int, int> new_pq;

        // TODO: this can probably be done more efficiently
        while (!importance_pq.empty()) {
            int contracted_importance = importance_pq.top_key();
            int contracted_node_id = importance_pq.top_value();
            importance_pq.pop();

            if (marked[contracted_node_id]) {
                new_pq.emplace(contracted_importance, contracted_node_id);
                continue;
            }

            if (cur_level != 0) {
                int new_importance = 0;
                switch (heuristic) {
                    case Heuristic::IN_OUT:
                        new_importance = inOutProductHeuristic(contracted, contracted_node_id);
                        break;
                    case Heuristic::EDGE_DIFFERENCE:
                        new_importance = edgeDifferenceHeuristic(contracted, contracted_node_id);
                        break;
                    case Heuristic::WEIGHTED_COST:
                        new_importance = weightedCostHeuristic(contracted, contracted_node_id);
                        break;
                    case Heuristic::MIXED:
                        new_importance = mixedHeuristic(contracted, contracted_node_id, cur_level);
                        break;
                    default:
                        break;
                }

                while (!importance_pq.empty() && new_importance > contracted_importance) {
                    if (marked[contracted_node_id])
                        new_pq.emplace(new_importance, contracted_node_id);
                    else
                        importance_pq.emplace(new_importance, contracted_node_id);
                    contracted_importance = importance_pq.top_key();
                    contracted_node_id = importance_pq.top_value();
                    importance_pq.pop();

                    switch (heuristic) {
                        case Heuristic::IN_OUT:
                            new_importance = inOutProductHeuristic(contracted, contracted_node_id);
                            break;
                        case Heuristic::EDGE_DIFFERENCE:
                            new_importance = edgeDifferenceHeuristic(contracted, contracted_node_id);
                            break;
                        case Heuristic::WEIGHTED_COST:
                            new_importance = weightedCostHeuristic(contracted, contracted_node_id);
                            break;
                        case Heuristic::MIXED:
                            new_importance = mixedHeuristic(contracted, contracted_node_id, cur_level);
                            break;
                        default:
                            break;
                    }
                }
            }

            if (marked[contracted_node_id])
                new_pq.emplace(contracted_importance, contracted_node_id);
            else {
                marked[contracted_node_id] = true;
                independent_set.push_back(contracted_node_id);
                for (ContractionEdge& e : m_graph_contr[contracted_node_id]) marked[e.m_target] = true;
                for (ContractionEdge& e : m_reverse_graph_contr[contracted_node_id]) marked[e.m_target] = true;
            }
        }
        importance_pq = std::move(new_pq);

// contract nodes from independent set in parallel
#pragma omp parallel for
        for (int i = 0; i < independent_set.size(); ++i) {
            contractNode(contracted, independent_set[i], omp_get_thread_num());
            contracted[independent_set[i]] = true;
            m_node_level[independent_set[i]] = cur_level;
        }

        for (int i = 0; i < m_num_threads; ++i) {
            for (auto& fwd_shortcut : m_contr_data[i].m_shortcuts_fwd)
                m_graph_contr[fwd_shortcut.first].push_back(fwd_shortcut.second);
            for (auto& bwd_shortcut : m_contr_data[i].m_shortcuts_bwd)
                m_reverse_graph_contr[bwd_shortcut.first].push_back(bwd_shortcut.second);

            m_contr_data[i].m_shortcuts_fwd.clear();
            m_contr_data[i].m_shortcuts_bwd.clear();

            if (i == 0) continue;

            for (int j = 0; j < m_num_nodes; ++j) {
                m_contr_data[0].m_num_contracted_neighbors[j] += m_contr_data[i].m_num_contracted_neighbors[j];
                m_contr_data[i].m_num_contracted_neighbors[j] = 0;
            }
        }

        num_contracted += static_cast<int>(independent_set.size());

        if (num_contracted >= m_num_nodes * 0.8) {
            omp_set_num_threads(1);
        }

        ++cur_level;
        std::cout << "Finished contracting: " << num_contracted << "\n";
    }

    m_contr_data.clear();
    omp_set_num_threads(m_num_threads);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    std::cout << "Finished creating CH. Took " << elapsed.count() << " seconds" << "\n";
}

void Graph::createContractionGraphs() {
    m_graph_contr.resize(m_num_nodes);
    m_reverse_graph_contr.resize(m_num_nodes);

    for (int i = 0; i < m_num_nodes; ++i) {
        uint64_t start_idx = m_graph_indices[i];
        uint64_t end_idx = m_graph_indices[i + 1];
        for (uint64_t j = start_idx; j < end_idx; ++j) {
            const Edge& e = m_graph_edges[j];
            m_graph_contr[i].push_back(ContractionEdge(e.m_target, e.m_cost));
            m_reverse_graph_contr[e.m_target].push_back(ContractionEdge(i, e.m_cost));
        }
    }
}

int Graph::inOutProductHeuristic(std::vector<bool>& contracted, int node) {
    int num_outgoing = 0;
    int num_incoming = 0;

    for (auto& outgoing : m_graph_contr[node]) {
        if (!contracted[outgoing.m_target]) ++num_outgoing;
    }

    for (auto& incoming : m_reverse_graph_contr[node]) {
        if (!contracted[incoming.m_target]) ++num_incoming;
    }
    return num_outgoing * num_incoming;
}

int Graph::edgeDifferenceHeuristic(std::vector<bool>& contracted, int node) {
    int num_edges_deleted = 0;
    for (auto& outgoing : m_graph_contr[node]) {
        if (!contracted[outgoing.m_target]) ++num_edges_deleted;
    }

    for (auto& incoming : m_reverse_graph_contr[node]) {
        if (!contracted[incoming.m_target]) ++num_edges_deleted;
    }

    int num_added_shortcuts = 0;

    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra
    int num_outgoing = 0;
    int max_distance_out = -1;
    for (auto& outgoing : m_graph_contr[node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        ++num_outgoing;
        m_contr_data[0].m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data[0].m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph_contr[node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance, 0);

        for (auto& outgoing : m_graph_contr[node]) {
            if (m_contr_data[0].m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data[0].m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data[0].m_visited[outgoing.m_target] = true;
            m_contr_data[0].m_reset_visited.push_back(outgoing.m_target);

            ++num_added_shortcuts;
        }

        // reset contraction data
        for (int& num : m_contr_data[0].m_reset_visited) m_contr_data[0].m_visited[num] = false;
        for (int& num : m_contr_data[0].m_reset_distances)
            m_contr_data[0].m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data[0].m_reset_visited.clear();
        m_contr_data[0].m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data[0].m_reset_outgoing) m_contr_data[0].m_outgoing[num] = false;
    m_contr_data[0].m_reset_outgoing.clear();

    return num_added_shortcuts - num_edges_deleted;
}

int Graph::weightedCostHeuristic(std::vector<bool>& contracted, int node) {
    int num_outgoing = 0;
    int num_incoming = 0;

    for (auto& outgoing : m_graph_contr[node]) {
        if (!contracted[outgoing.m_target]) ++num_outgoing;
    }

    for (auto& incoming : m_reverse_graph_contr[node]) {
        if (!contracted[incoming.m_target]) ++num_incoming;
    }

    int max_cost = 0;
    int num_added_shortcuts = 0;

    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra

    int max_distance_out = -1;
    for (auto& outgoing : m_graph_contr[node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        m_contr_data.at(0).m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data[0].m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph_contr[node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance, 0);

        for (auto& outgoing : m_graph_contr[node]) {
            if (m_contr_data[0].m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data[0].m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data[0].m_visited[outgoing.m_target] = true;
            m_contr_data[0].m_reset_visited.push_back(outgoing.m_target);

            ++num_added_shortcuts;

            if (incoming.m_cost + outgoing.m_cost > max_cost) max_cost = incoming.m_cost + outgoing.m_cost;
        }

        // reset contraction data
        for (int& num : m_contr_data[0].m_reset_visited) m_contr_data[0].m_visited[num] = false;
        for (int& num : m_contr_data[0].m_reset_distances)
            m_contr_data[0].m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data[0].m_reset_visited.clear();
        m_contr_data[0].m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data[0].m_reset_outgoing) m_contr_data[0].m_outgoing[num] = false;
    m_contr_data[0].m_reset_outgoing.clear();

    return 0.8 * max_cost + 0.2 * num_outgoing * num_incoming;
    // return max_cost;
}

int Graph::mixedHeuristic(std::vector<bool>& contracted, int node, int cur_level) {
    int num_outgoing = 0;
    int num_incoming = 0;

    for (auto& outgoing : m_graph_contr[node]) {
        if (!contracted[outgoing.m_target]) ++num_outgoing;
    }

    for (auto& incoming : m_reverse_graph_contr[node]) {
        if (!contracted[incoming.m_target]) ++num_incoming;
    }

    int max_cost = 0;
    int num_added_shortcuts = 0;

    int max_neighbour_level = -1;

    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra
    int max_distance_out = -1;
    for (auto& outgoing : m_graph_contr[node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        m_contr_data[0].m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data[0].m_outgoing[outgoing.m_target] = true;
    }

    int underlying_shortcuts = 0;

    for (auto& incoming : m_reverse_graph_contr[node]) {
        if (contracted[incoming.m_target]) {
            if (m_node_level[incoming.m_target] > max_neighbour_level)
                max_neighbour_level = m_node_level[incoming.m_target];
            continue;
        }
        int max_distance = incoming.m_cost + max_distance_out;
        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance, 0);

        for (auto& outgoing : m_graph_contr[node]) {
            if (m_contr_data[0].m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) {
                if (m_node_level[incoming.m_target] > max_neighbour_level)
                    max_neighbour_level = m_node_level[incoming.m_target];
                continue;
            }
            if (m_contr_data[0].m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data[0].m_visited[outgoing.m_target] = true;
            m_contr_data[0].m_reset_visited.push_back(outgoing.m_target);

            if (incoming.m_num_underlying_arcs != -1) underlying_shortcuts += incoming.m_num_underlying_arcs;
            if (outgoing.m_num_underlying_arcs != -1) underlying_shortcuts += outgoing.m_num_underlying_arcs;

            ++num_added_shortcuts;

            if (incoming.m_cost + outgoing.m_cost > max_cost) max_cost = incoming.m_cost + outgoing.m_cost;
        }

        // reset contraction data
        for (int& num : m_contr_data[0].m_reset_visited) m_contr_data[0].m_visited[num] = false;
        for (int& num : m_contr_data[0].m_reset_distances)
            m_contr_data[0].m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data[0].m_reset_visited.clear();
        m_contr_data[0].m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data[0].m_reset_outgoing) m_contr_data[0].m_outgoing[num] = false;
    m_contr_data[0].m_reset_outgoing.clear();

    if (max_neighbour_level > cur_level)
        max_neighbour_level = 0;
    else
        ++max_neighbour_level;

    return static_cast<int>(0.001 * static_cast<double>(max_cost)) +
           2 * (num_added_shortcuts - (num_incoming + num_outgoing)) +
           1 * m_contr_data[0].m_num_contracted_neighbors[node] + 5 * max_neighbour_level + underlying_shortcuts;
}

void Graph::contractNode(std::vector<bool>& contracted, int contracted_node, int thread_num) {
    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra
    int num_outgoing = 0;
    int max_distance_out = -1;

    for (auto& outgoing : m_graph_contr[contracted_node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        ++m_contr_data[thread_num].m_num_contracted_neighbors[outgoing.m_target];

        ++num_outgoing;
        m_contr_data[thread_num].m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data[thread_num].m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph_contr[contracted_node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        ++m_contr_data[thread_num].m_num_contracted_neighbors[incoming.m_target];

        contractionDijkstra(incoming.m_target, contracted_node, contracted, num_outgoing, max_distance, thread_num);

        for (auto& outgoing : m_graph_contr[contracted_node]) {
            if (m_contr_data[thread_num].m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data[thread_num].m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data[thread_num].m_visited[outgoing.m_target] = true;
            m_contr_data[thread_num].m_reset_visited.push_back(outgoing.m_target);

            // add shortcut
            int underlying_arcs = 0;
            if (incoming.m_num_underlying_arcs != -1) underlying_arcs += incoming.m_num_underlying_arcs;
            if (outgoing.m_num_underlying_arcs != -1) underlying_arcs += outgoing.m_num_underlying_arcs;

            ContractionEdge fwd_shortcut(outgoing.m_target, incoming.m_cost + outgoing.m_cost, contracted_node,
                                         underlying_arcs);
            m_contr_data[thread_num].m_shortcuts_fwd.push_back(std::make_pair(incoming.m_target, fwd_shortcut));
            ContractionEdge bwd_shortcut(incoming.m_target, incoming.m_cost + outgoing.m_cost, contracted_node,
                                         underlying_arcs);
            m_contr_data[thread_num].m_shortcuts_bwd.push_back(std::make_pair(outgoing.m_target, bwd_shortcut));
        }

        // reset contraction data
        for (int& num : m_contr_data[thread_num].m_reset_visited) m_contr_data[thread_num].m_visited[num] = false;
        for (int& num : m_contr_data[thread_num].m_reset_distances)
            m_contr_data[thread_num].m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data[thread_num].m_reset_visited.clear();
        m_contr_data[thread_num].m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data[thread_num].m_reset_outgoing) m_contr_data[thread_num].m_outgoing[num] = false;
    m_contr_data[thread_num].m_reset_outgoing.clear();
}

void Graph::contractionDijkstra(int start, int contracted_node, std::vector<bool>& contracted, int num_outgoing,
                                int max_distance, int thread_num) {
    radix_heap::pair_radix_heap<int, int> pq;

    int num_visited_outgoing = 0;
    m_contr_data[thread_num].m_distances[start] = 0;
    m_contr_data[thread_num].m_reset_distances.push_back(start);
    pq.emplace(0, start);

    while (!pq.empty()) {
        int cur_dist = pq.top_key();
        int cur_node_id = pq.top_value();
        pq.pop();

        if (m_contr_data[thread_num].m_distances[cur_node_id] != cur_dist) continue;

        if (m_contr_data[thread_num].m_outgoing[cur_node_id]) ++num_visited_outgoing;
        if (cur_dist > max_distance || num_visited_outgoing == num_outgoing) break;

        for (ContractionEdge& e : m_graph_contr[cur_node_id]) {
            if (contracted[e.m_target]) continue;
            if (m_contr_data[thread_num].m_distances[e.m_target] >
                m_contr_data[thread_num].m_distances[cur_node_id] + e.m_cost) {
                if (m_contr_data[thread_num].m_distances[e.m_target] == std::numeric_limits<int>::max())
                    m_contr_data[thread_num].m_reset_distances.push_back(e.m_target);
                m_contr_data[thread_num].m_distances[e.m_target] =
                    m_contr_data[thread_num].m_distances[cur_node_id] + e.m_cost;
                pq.emplace(m_contr_data[thread_num].m_distances[e.m_target], e.m_target);
            }
        }
    }
}

int Graph::simplifiedHubLabelQuery(std::vector<std::tuple<uint32_t, uint32_t, uint16_t, uint16_t>>& fwd_labels,
                                   int node) {
    int distance = std::numeric_limits<int>::max();
    auto fwd_iter = fwd_labels.begin();
    uint64_t bwd_node_index = m_bwd_indices[node];
    uint64_t bwd_next_index = m_bwd_indices[node + 1];

    while (fwd_iter != fwd_labels.end() && bwd_node_index < bwd_next_index) {
        if (std::get<0>(*fwd_iter) == std::get<0>(m_bwd_hub_labels[bwd_node_index])) {
            int total_dist = static_cast<int>(std::get<1>(*fwd_iter)) +
                             static_cast<int>(std::get<1>(m_bwd_hub_labels[bwd_node_index]));
            if (total_dist < distance) distance = total_dist;
            ++fwd_iter;
            ++bwd_node_index;
        } else if (std::get<0>(*fwd_iter) < std::get<0>(m_bwd_hub_labels[bwd_node_index])) {
            ++fwd_iter;
        } else {
            ++bwd_node_index;
        }
    }

    return distance;
}

int Graph::simplifiedHubLabelQuery(int node,
                                   std::vector<std::tuple<uint32_t, uint32_t, uint16_t, uint16_t>>& bwd_labels) {
    int distance = std::numeric_limits<int>::max();
    uint64_t fwd_node_index = m_fwd_indices[node];
    uint64_t fwd_next_index = m_fwd_indices[node + 1];
    auto bwd_iter = bwd_labels.begin();

    while (fwd_node_index < fwd_next_index && bwd_iter != bwd_labels.end()) {
        if (std::get<0>(*bwd_iter) == std::get<0>(m_fwd_hub_labels[fwd_node_index])) {
            int total_dist = static_cast<int>(std::get<1>(*bwd_iter)) +
                             static_cast<int>(std::get<1>(m_fwd_hub_labels[fwd_node_index]));
            if (total_dist < distance) distance = total_dist;
            ++bwd_iter;
            ++fwd_node_index;
        } else if (std::get<0>(*bwd_iter) < std::get<0>(m_fwd_hub_labels[fwd_node_index])) {
            ++bwd_iter;
        } else {
            ++fwd_node_index;
        }
    }

    return distance;
}

void Graph::readGraphFromCHFMI(const std::string& path) {
    std::cout << "Started reading graph file from CHFMI." << "\n";

    auto begin = std::chrono::high_resolution_clock::now();
    std::ifstream infile;
    try {
        infile.open(path);
        if (!infile.good()) throw std::runtime_error("File doesn\\'t exist!");
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        exit(1);
    }

    std::string line = "#";

    // skip the metadata which begins with #
    while (line[0] == '#') getline(infile, line);

    getline(infile, line);
    m_num_nodes = std::stoi(line);
    getline(infile, line);
    int num_edges = std::stoi(line);

    m_node_level.resize(m_num_nodes);
    m_node_coords.resize(m_num_nodes);
    m_graph_contr.resize(m_num_nodes);

    for (int i = 0; i < m_num_nodes; ++i) {
        getline(infile, line);
        std::stringstream ss(line);
        std::string s;
        int node_id;
        int ch_id;
        double lat, lon;

        getline(ss, s, ' ');
        node_id = std::stoi(s);
        getline(ss, s, ' ');
        ch_id = std::stoi(s);
        getline(ss, s, ' ');
        lon = stod(s);
        getline(ss, s, ' ');
        lat = stod(s);

        m_node_coords[node_id] = {lon, lat};
        m_node_level[node_id] = ch_id;
    }

    int edge_count = 0;

    for (int i = 0; i < num_edges; ++i) {
        getline(infile, line);
        std::stringstream ss(line);
        std::string s;
        int src_id, tgt_id, cost, contraction_node_id;

        try {
            getline(ss, s, ' ');
            src_id = std::stoi(s);
            getline(ss, s, ' ');
            tgt_id = std::stoi(s);
            getline(ss, s, ' ');
            cost = std::stoi(s);
            getline(ss, s, ' ');
            contraction_node_id = std::stoi(s);
            m_graph_contr[src_id].push_back(ContractionEdge(tgt_id, cost, contraction_node_id, -1));
            edge_count++;
        } catch (const std::exception& e) {
            std::cerr << "Error parsing edge line: " << line << "\n";
            std::cerr << "Exception: " << e.what() << "\n";
            continue;
        }
    }

    std::cout << "Edge Count: " << edge_count << "\n";

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished reading graph file from CHFMI. Took " << elapsed.count() << " milliseconds " << "\n";
}

void Graph::writeGraphToCHFMI(const std::string& filename) {
    std::ofstream out(filename);
    out << "# Timestamp: " << std::time(nullptr) << '\n';
    out << "# Type: Coastlines with CH Information\n\n";

    int num_edges = 0;
    for (const auto& edges : m_graph_contr) {
        num_edges += edges.size();
    }

    out << m_num_nodes << '\n';
    out << num_edges << '\n';

    for (int i = 0; i < m_num_nodes; ++i) {
        const auto& [lon, lat] = m_node_coords[i];
        // Assuming m_node_level contains the CH ID (rank)
        out << i << " " << m_node_level[i] << " " << lon << " " << lat << " 0" << '\n';
    }

    for (int i = 0; i < m_num_nodes; ++i) {
        const auto& edges = m_graph_contr[i];
        for (const auto& edge : edges) {
            out << i << " " << edge.m_target << " " << edge.m_cost << " " << edge.m_contraction_node << " 0" << '\n';
        }
    }
    out.close();
}

}  // namespace labosm
