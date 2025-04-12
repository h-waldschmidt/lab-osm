#include "graph.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <unordered_set>

#define _USE_MATH_DEFINES
#include <math.h>

namespace labosm {

Graph::Graph(const std::string& path, bool ch_available, Heuristic ch_heuristic)
    : m_ch_available(ch_available), m_is(false), m_num_threads(1) {
    readGraph(path);
    m_num_nodes = m_graph.size();

    if (m_ch_available) {
        createCHwithoutIS(ch_heuristic);
    }

    createReverseGraph();
}

Graph::Graph(const std::string& path, bool ch_available, int num_threads, Heuristic ch_heuristic)
    : m_ch_available(ch_available), m_num_threads(num_threads), m_is(true) {
    readGraph(path);
    m_num_nodes = m_graph.size();

    omp_set_num_threads(m_num_threads);

    if (m_ch_available) {
        createCHwithIS(ch_heuristic);
    }

    createReverseGraph();
}

void Graph::readGraph(const std::string& path) {
    std::cout << "Started reading graph file." << std::endl;

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

    // save node rank if CH Mode, else skip the node information
    // or save coords if distances in meters are required

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

    m_graph.clear();
    m_graph.resize(num_nodes);

    // read edge information
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

        m_graph[src].push_back(Edge{target, cost});
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished reading graph file. Took " << elapsed.count() << " milliseconds " << std::endl;
}

void Graph::createReverseGraph() {
    std::cout << "Started creating reverse graph." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty()) {
        std::cout << "Can't create reverse graph, because graph is empty" << std::endl;
        return;
    }

    m_reverse_graph.clear();
    m_reverse_graph.resize(m_num_nodes);

    if (m_ch_available)
        createReverseGraphCH();
    else
        createReverseGraphNormal();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating reverse graph. Took " << elapsed.count() << " milliseconds " << std::endl;
}

void Graph::dijkstraQuery(DijkstraQueryData& data) {
    if (data.needReset()) {
        data.reset();
    }

    struct PQEntry {
        int m_target;
        int m_distance;
        int m_previous;

        bool operator>(const PQEntry& other) const { return m_distance > other.m_distance; }
    };

    std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

    pq.push({data.m_start, 0, data.m_start});

    PQEntry element;
    while (!pq.empty()) {
        element = pq.top();
        pq.pop();
        data.num_pq_pops++;

        if (data.m_distances[element.m_target] <= element.m_distance) continue;

        data.m_distances[element.m_target] = element.m_distance;
        data.m_prev[element.m_target] = element.m_previous;
        data.m_reset_nodes.push_back(element.m_target);

        if (element.m_target == data.m_end) {
            data.m_distance = element.m_distance;
            return;
        }

        for (const Edge& edge : m_graph[element.m_target]) {
            int new_distance = element.m_distance + edge.m_cost;
            if (new_distance < data.m_distances[edge.m_target]) {
                pq.push(PQEntry{edge.m_target, new_distance, element.m_target});
            }
        }
    }
    data.m_distance = data.m_distances[data.m_end];

    // TODO: return path
}

void Graph::bidirectionalDijkstraGetPath(QueryData& data) {
    if (!(data.m_distance < std::numeric_limits<int>::max() || data.m_distance > -1)) {
        std::cout << "Can't return path for invalid data!" << std::endl;
        return;
    }

    data.m_shortest_path.clear();

    std::vector<int> fwd_path;
    std::vector<int> bwd_path;

    int cur_node = data.m_meeting_node;
    while (cur_node != data.m_start) {
        cur_node = data.m_fwd_prev[cur_node];
        fwd_path.push_back(cur_node);
    }

    cur_node = data.m_meeting_node;
    while (cur_node != data.m_end) {
        cur_node = data.m_bwd_prev[cur_node];
        bwd_path.push_back(cur_node);
    }

    for (int i = fwd_path.size() - 1; i >= 0; --i) {
        data.m_shortest_path.push_back(fwd_path[i]);
    }

    data.m_shortest_path.push_back(data.m_meeting_node);

    for (int i = 0; i < bwd_path.size(); ++i) {
        data.m_shortest_path.push_back(bwd_path[i]);
    }
    std::cout << "Level of Meeting Node: " << m_node_level[data.m_meeting_node] << std::endl;
    for (size_t i = 0; i < data.m_shortest_path.size(); i++) {
        std::cout << "NodeID: " << data.m_shortest_path[i] << " Level: " << m_node_level[data.m_shortest_path[i]]
                  << std::endl;
    }
}

void Graph::contractionHierarchyQuery(QueryData& data) {
    if (data.m_start < 0 || data.m_end < 0) {
        std::cout << "Invalid start or end nodes!" << std::endl;
        return;
    }

    if (data.needReset()) {
        data.reset();
    }

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        fwd_pq;
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        bwd_pq;

    data.m_distances_fwd[data.m_start] = 0;
    data.m_distances_bwd[data.m_end] = 0;

    // first corresponds to distance and second is node index
    fwd_pq.push(std::make_pair(0, data.m_start));
    bwd_pq.push(std::make_pair(0, data.m_end));

    data.m_visited_fwd[data.m_start] = true;
    data.m_visited_bwd[data.m_end] = true;

    data.m_reset_nodes_fwd.push_back(data.m_start);
    data.m_reset_nodes_bwd.push_back(data.m_end);

    std::pair<int, int> fwd_node;
    std::pair<int, int> bwd_node;

    while (!fwd_pq.empty() || !bwd_pq.empty()) {
        while (!fwd_pq.empty()) {
            fwd_node = fwd_pq.top();
            fwd_pq.pop();
            data.num_pq_pops++;

            if (data.m_visited_fwd[fwd_node.second] && fwd_node.first > data.m_distances_fwd[fwd_node.second]) continue;
            if (fwd_node.first > data.m_distance) break;

            // forward step
            for (Edge& e : m_graph[fwd_node.second]) {
                if (!data.m_visited_fwd[e.m_target] || data.m_distances_fwd[e.m_target] > fwd_node.first + e.m_cost) {
                    if (!data.m_visited_fwd[e.m_target]) data.m_reset_nodes_fwd.push_back(e.m_target);

                    data.m_distances_fwd[e.m_target] = fwd_node.first + e.m_cost;
                    fwd_pq.push(std::make_pair(data.m_distances_fwd[e.m_target], e.m_target));
                    data.m_visited_fwd[e.m_target] = true;

                    data.m_fwd_prev[e.m_target] = fwd_node.second;
                }

                if (data.m_visited_bwd[e.m_target] &&
                    data.m_distances_fwd[fwd_node.second] + e.m_cost + data.m_distances_bwd[e.m_target] <
                        data.m_distance) {
                    data.m_distance =
                        data.m_distances_fwd[fwd_node.second] + e.m_cost + data.m_distances_bwd[e.m_target];
                    data.m_meeting_node = e.m_target;
                }
            }

            break;
        }

        while (!bwd_pq.empty()) {
            bwd_node = bwd_pq.top();
            bwd_pq.pop();
            data.num_pq_pops++;

            if (data.m_visited_bwd[bwd_node.second] && bwd_node.first > data.m_distances_fwd[bwd_node.second]) continue;
            if (bwd_node.first > data.m_distance) break;

            // backward step
            for (Edge& e : m_reverse_graph[bwd_node.second]) {
                if (!data.m_visited_bwd[e.m_target] || data.m_distances_bwd[e.m_target] > bwd_node.first + e.m_cost) {
                    if (!data.m_visited_fwd[e.m_target]) data.m_reset_nodes_bwd.push_back(e.m_target);

                    data.m_distances_bwd[e.m_target] = bwd_node.first + e.m_cost;
                    bwd_pq.push(std::make_pair(data.m_distances_bwd[e.m_target], e.m_target));
                    data.m_visited_bwd[e.m_target] = true;

                    data.m_bwd_prev[e.m_target] = bwd_node.second;
                }

                if (data.m_visited_fwd[e.m_target] &&
                    data.m_distances_bwd[bwd_node.second] + e.m_cost + data.m_distances_fwd[e.m_target] <
                        data.m_distance) {
                    data.m_distance =
                        data.m_distances_bwd[bwd_node.second] + e.m_cost + data.m_distances_fwd[e.m_target];
                    data.m_meeting_node = e.m_target;
                }
            }
            break;
        }

        // termination condition
        if (fwd_node.first >= data.m_distance || bwd_node.first >= data.m_distance) {
            break;
        }
    }
}

void Graph::createHubLabelsWithoutIS(int threshold) {
    std::cout << "Started creating hub labels." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty() || m_reverse_graph.empty()) {
        std::cout << "Can't create hub labels, because graph or reverse graph is empty" << std::endl;
        return;
    }

    m_fwd_indices.clear();
    m_fwd_indices.resize(m_num_nodes + 1, 0);
    m_bwd_indices.clear();
    m_bwd_indices.resize(m_num_nodes + 1, 0);

    m_fwd_hub_labels.clear();
    m_bwd_hub_labels.clear();

    // sort the levels, but don't change the original vector
    m_level_indices_sorted.clear();
    m_level_indices_sorted.resize(m_num_nodes);
    std::iota(m_level_indices_sorted.begin(), m_level_indices_sorted.end(), 0);
    std::sort(m_level_indices_sorted.begin(), m_level_indices_sorted.end(),
              [&](int i, int j) { return m_node_level[i] > m_node_level[j]; });

    m_node_indices.clear();
    m_node_indices.resize(m_num_nodes);
    // save for each node its corresponding index
    for (int i = 0; i < m_num_nodes; ++i) {
        m_node_indices[m_level_indices_sorted[i]] = i;
    }

    std::vector<std::pair<int, int>> fwd_labels;
    std::vector<std::pair<int, int>> bwd_labels;

    for (int i = 0; i < m_num_nodes; ++i) {
        int node = m_level_indices_sorted[i];

        if (i == 0) {
            m_fwd_hub_labels.push_back(std::make_pair(node, 0));
            m_bwd_hub_labels.push_back(std::make_pair(node, 0));
            m_fwd_indices[1] = 1;
            m_bwd_indices[1] = 1;
            continue;
        }

        fwd_labels.clear();
        bwd_labels.clear();

        fwd_labels.push_back(std::make_pair(node, 0));
        bwd_labels.push_back(std::make_pair(node, 0));

        int test_1 = m_node_indices[node];
        int test_2 = m_level_indices_sorted[i];

        // fwd labels
        for (Edge& e : m_graph[node]) {
            if (m_node_level[node] >= m_node_level[e.m_target]) continue;

            for (uint32_t j = m_fwd_indices[m_node_indices[e.m_target]];
                 j < m_fwd_indices[m_node_indices[e.m_target] + 1]; j++) {
                if (m_fwd_hub_labels[j].second + e.m_cost <= threshold)
                    fwd_labels.push_back(
                        std::make_pair(m_fwd_hub_labels[j].first, m_fwd_hub_labels[j].second + e.m_cost));
            }
        }
        // remove duplicates
        std::sort(fwd_labels.begin(), fwd_labels.end(),
                  [](auto& left, auto& right) { return left.first < right.first; });
        for (auto iter = fwd_labels.begin(); iter != fwd_labels.end();) {
            auto iter_2 = fwd_labels.end() - 1;
            if (std::distance(iter, iter_2) != 0 && iter->first == (iter + 1)->first) {
                if (iter->second >= (iter + 1)->second) {
                    iter = fwd_labels.erase(iter);
                } else {
                    iter = fwd_labels.erase(iter + 1);
                    --iter;
                }
            } else {
                ++iter;
            }
        }

        for (auto iter = fwd_labels.begin(); iter != fwd_labels.end();) {
            int best_dist = std::numeric_limits<int>::max();
            if (iter->first != node) best_dist = simplifiedHubLabelQuery(fwd_labels, iter->first);
            if (best_dist < iter->second)
                iter = fwd_labels.erase(iter);
            else
                ++iter;
        }

        // bwd labels
        for (Edge& e : m_reverse_graph[node]) {
            if (m_node_level[node] >= m_node_level[e.m_target]) continue;

            for (uint32_t j = m_bwd_indices[m_node_indices[e.m_target]];
                 j < m_bwd_indices[m_node_indices[e.m_target] + 1]; j++) {
                if (m_bwd_hub_labels[j].second + e.m_cost <= threshold)
                    bwd_labels.push_back(
                        std::make_pair(m_bwd_hub_labels[j].first, m_bwd_hub_labels[j].second + e.m_cost));
            }
        }
        // remove duplicates
        std::sort(bwd_labels.begin(), bwd_labels.end(),
                  [](auto& left, auto& right) { return left.first < right.first; });
        for (auto iter = bwd_labels.begin(); iter != bwd_labels.end();) {
            auto iter_2 = bwd_labels.end() - 1;
            if (std::distance(iter, iter_2) != 0 && iter->first == (iter + 1)->first) {
                if (iter->second >= (iter + 1)->second) {
                    iter = bwd_labels.erase(iter);
                } else {
                    iter = bwd_labels.erase(iter + 1);
                    --iter;
                }
            } else {
                ++iter;
            }
        }
        for (auto iter = bwd_labels.begin(); iter != bwd_labels.end();) {
            int best_dist = std::numeric_limits<int>::max();
            if (iter->first != node) best_dist = simplifiedHubLabelQuery(iter->first, bwd_labels);
            if (best_dist < iter->second)
                iter = bwd_labels.erase(iter);
            else
                ++iter;
        }
        // update hub label data
        m_fwd_indices[i + 1] = m_fwd_indices[i] + static_cast<uint32_t>(fwd_labels.size());
        m_bwd_indices[i + 1] = m_bwd_indices[i] + static_cast<uint32_t>(bwd_labels.size());

        for (auto label : fwd_labels) {
            m_fwd_hub_labels.push_back(label);
        }
        for (auto label : bwd_labels) {
            m_bwd_hub_labels.push_back(label);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating hub labels. Took " << elapsed.count() << " milliseconds " << std::endl;
}

void Graph::createHubLabelsWithIS(int threshold) {
    std::cout << "Started creating hub labels." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty() || m_reverse_graph.empty()) {
        std::cout << "Can't create hub labels, because graph or reverse graph is empty" << std::endl;
        return;
    }

    m_fwd_indices.clear();
    m_fwd_indices.resize(m_num_nodes + 1, 0);
    m_bwd_indices.clear();
    m_bwd_indices.resize(m_num_nodes + 1, 0);

    m_fwd_hub_labels.clear();
    m_bwd_hub_labels.clear();

    // sort the levels, but don't change the original vector
    std::unordered_set<int> different_levels;
    for (int i = 0; i < m_num_nodes; i++) {
        different_levels.emplace(m_node_level[i]);
    }
    std::vector<int> different_levels_vec;
    different_levels_vec.insert(different_levels_vec.end(), different_levels.begin(), different_levels.end());
    std::sort(different_levels_vec.begin(), different_levels_vec.end(),
              [](int& left, int& right) { return left > right; });
    std::vector<std::vector<int>> level_buckets(different_levels_vec.size());
    for (int i = 0; i < m_num_nodes; ++i) {
        auto iter = find(different_levels_vec.begin(), different_levels_vec.end(), m_node_level[i]);
        level_buckets[iter - different_levels_vec.begin()].push_back(i);
    }

    // sort the levels, but don't change the original vector
    m_level_indices_sorted.clear();

    for (int i = 0; i < different_levels_vec.size(); ++i) {
        for (int j = 0; j < level_buckets[i].size(); ++j) {
            m_level_indices_sorted.push_back(level_buckets[i][j]);
        }
    }

    m_node_indices.clear();
    m_node_indices.resize(m_num_nodes);
    // save for each node its corresponding index
    for (int i = 0; i < m_num_nodes; ++i) {
        m_node_indices[m_level_indices_sorted[i]] = i;
    }

    uint32_t num_calculated = 0;

    for (int i = 0; i < different_levels_vec.size(); ++i) {
        if (i == 0) {
            int node = level_buckets[i][0];
            m_fwd_hub_labels.push_back(std::make_pair(node, 0));
            m_bwd_hub_labels.push_back(std::make_pair(node, 0));
            m_fwd_indices[1] = 1;
            m_bwd_indices[1] = 1;
            continue;
        }

#pragma omp parallel for ordered schedule(static, 1)
        for (int j = 0; j < level_buckets[i].size(); ++j) {
            int node = level_buckets[i][j];

            std::vector<std::pair<int, int>> fwd_labels;
            std::vector<std::pair<int, int>> bwd_labels;

            fwd_labels.push_back(std::make_pair(node, 0));
            bwd_labels.push_back(std::make_pair(node, 0));

            // fwd labels
            for (Edge& e : m_graph[node]) {
                if (m_node_level[node] >= m_node_level[e.m_target]) continue;

                for (uint32_t j = m_fwd_indices[m_node_indices[e.m_target]];
                     j < m_fwd_indices[m_node_indices[e.m_target] + 1]; j++) {
                    if (m_fwd_hub_labels[j].second + e.m_cost <= threshold)
                        fwd_labels.push_back(
                            std::make_pair(m_fwd_hub_labels[j].first, m_fwd_hub_labels[j].second + e.m_cost));
                }
            }

            // remove duplicates
            std::sort(fwd_labels.begin(), fwd_labels.end(),
                      [](auto& left, auto& right) { return left.first < right.first; });
            for (auto iter = fwd_labels.begin(); iter != fwd_labels.end();) {
                auto iter_2 = fwd_labels.end() - 1;
                if (std::distance(iter, iter_2) != 0 && iter->first == (iter + 1)->first) {
                    if (iter->second >= (iter + 1)->second) {
                        iter = fwd_labels.erase(iter);
                    } else {
                        iter = fwd_labels.erase(iter + 1);
                        --iter;
                    }
                } else {
                    ++iter;
                }
            }

            for (auto iter = fwd_labels.begin(); iter != fwd_labels.end();) {
                int best_dist = std::numeric_limits<int>::max();
                if (iter->first != node) best_dist = simplifiedHubLabelQuery(fwd_labels, iter->first);
                if (best_dist < iter->second)
                    iter = fwd_labels.erase(iter);
                else
                    ++iter;
            }

            // bwd labels
            for (Edge& e : m_reverse_graph[node]) {
                if (m_node_level[node] >= m_node_level[e.m_target]) continue;

                for (uint32_t j = m_bwd_indices[m_node_indices[e.m_target]];
                     j < m_bwd_indices[m_node_indices[e.m_target] + 1]; j++) {
                    if (m_bwd_hub_labels[j].second + e.m_cost <= threshold)
                        bwd_labels.push_back(
                            std::make_pair(m_bwd_hub_labels[j].first, m_bwd_hub_labels[j].second + e.m_cost));
                }
            }

            // remove duplicates
            std::sort(bwd_labels.begin(), bwd_labels.end(),
                      [](auto& left, auto& right) { return left.first < right.first; });
            for (auto iter = bwd_labels.begin(); iter != bwd_labels.end();) {
                auto iter_2 = bwd_labels.end() - 1;
                if (std::distance(iter, iter_2) != 0 && iter->first == (iter + 1)->first) {
                    if (iter->second >= (iter + 1)->second) {
                        iter = bwd_labels.erase(iter);
                    } else {
                        iter = bwd_labels.erase(iter + 1);
                        --iter;
                    }
                } else {
                    ++iter;
                }
            }

            for (auto iter = bwd_labels.begin(); iter != bwd_labels.end();) {
                int best_dist = std::numeric_limits<int>::max();
                if (iter->first != node) best_dist = simplifiedHubLabelQuery(iter->first, bwd_labels);
                if (best_dist < iter->second)
                    iter = bwd_labels.erase(iter);
                else
                    ++iter;
            }

#pragma omp ordered
            {
                int fwd_node_index = m_node_indices[node];
                m_fwd_indices[fwd_node_index + 1] =
                    m_fwd_indices[fwd_node_index] + static_cast<uint32_t>(fwd_labels.size());
                m_bwd_indices[fwd_node_index + 1] =
                    m_bwd_indices[fwd_node_index] + static_cast<uint32_t>(bwd_labels.size());

                for (auto label : fwd_labels) {
                    m_fwd_hub_labels.push_back(label);
                }
                for (auto label : bwd_labels) {
                    m_bwd_hub_labels.push_back(label);
                }
            }
        }
        num_calculated += static_cast<uint32_t>(level_buckets[i].size());
        // std::cout << "Finished Hub Labels for " << num_calculated << " num of nodes." << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating hub labels. Took " << elapsed.count() << " milliseconds " << std::endl;
}

void Graph::hubLabelQuery(QueryData& data) {
    if (data.m_start < 0 || data.m_end < 0) {
        std::cout << "Invalid start or end nodes!" << std::endl;
        return;
    }

    std::cout << "Started Hub Query." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    data.m_distance = std::numeric_limits<int>::max();
    data.m_meeting_node = -1;

    u_int64_t fwd_node_index = m_fwd_indices[m_node_indices[data.m_start]];
    u_int64_t fwd_next_index = m_fwd_indices[m_node_indices[data.m_start] + 1];

    u_int64_t bwd_node_index = m_bwd_indices[m_node_indices[data.m_end]];
    u_int64_t bwd_next_index = m_bwd_indices[m_node_indices[data.m_end] + 1];

    while (fwd_node_index < fwd_next_index && bwd_node_index < bwd_next_index) {
        if (m_fwd_hub_labels[fwd_node_index].first == m_bwd_hub_labels[bwd_node_index].first) {
            if (m_fwd_hub_labels[fwd_node_index].second + m_bwd_hub_labels[bwd_node_index].second < data.m_distance) {
                data.m_meeting_node = m_fwd_hub_labels[fwd_node_index].first;
                data.m_distance = m_fwd_hub_labels[fwd_node_index].second + m_bwd_hub_labels[bwd_node_index].second;
            }
            ++fwd_node_index;
            ++bwd_node_index;
        } else if (m_fwd_hub_labels[fwd_node_index].first < m_bwd_hub_labels[bwd_node_index].first) {
            ++fwd_node_index;
        } else {
            ++bwd_node_index;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "Finished Hub Label query. Took " << elapsed.count() << " nanoseconds " << std::endl;
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

std::vector<int> intersection(std::vector<int> v1, std::vector<int> v2) {
    std::vector<int> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));
    return v3;
}

void Graph::createReverseGraphNormal() {
    for (int i = 0; i < m_num_nodes; ++i) {
        for (Edge e : m_graph[i]) {
            int new_source = e.m_target;
            e.m_target = i;

            m_reverse_graph[new_source].push_back(e);
        }
    }
}

void Graph::createReverseGraphCH() {
    m_graph.clear();
    m_graph.resize(m_num_nodes);
    m_reverse_graph.clear();
    m_reverse_graph.resize(m_num_nodes);
    m_reverse_graph_contr.clear();
    m_reverse_graph_contr.resize(m_num_nodes);

    for (int i = 0; i < m_num_nodes; ++i) {
        for (auto iter = m_graph_contr[i].begin(); iter != m_graph_contr[i].end();) {
            if (m_node_level[i] > m_node_level[iter->m_target]) {
                ContractionEdge e = *iter;
                int new_source = e.m_target;
                e.m_target = i;
                m_reverse_graph_contr[new_source].push_back(e);
                m_graph_contr[i].erase(iter);
            } else {
                iter++;
            }
        }
    }

    // copy the data
    for (int i = 0; i < m_num_nodes; ++i) {
        for (const auto& e : m_graph_contr[i]) {
            m_graph[i].push_back(Edge(e.m_target, e.m_cost));
        }
        for (const auto& e : m_reverse_graph_contr[i]) {
            m_reverse_graph[i].push_back(Edge(e.m_target, e.m_cost));
        }
    }

    // find the correct shortcut pointers
    for (int i = 0; i < m_num_nodes; ++i) {
        for (int j = 0; j < m_graph_contr[i].size(); ++j) {
            const ContractionEdge& e = m_graph_contr[i][j];
            if (e.isShortcut()) {
                Edge* child_1 = nullptr;
                for (int k = 0; k < m_graph[i].size(); ++k) {
                    if (m_graph[i][k].m_target == e.m_contraction_node) {
                        child_1 = &m_graph[i][k];
                        break;
                    }
                }
                Edge* child_2 = nullptr;
                for (int k = 0; k < m_graph[e.m_contraction_node].size(); ++k) {
                    if (m_graph[e.m_contraction_node][k].m_target == e.m_target) {
                        child_2 = &m_graph[e.m_contraction_node][k];
                        break;
                    }
                }
                m_graph[i][j].m_child_1 = child_1;
                m_graph[i][j].m_child_2 = child_2;
            }
        }

        for (int j = 0; j < m_reverse_graph_contr[i].size(); ++j) {
            const ContractionEdge& e = m_reverse_graph_contr[i][j];
            if (e.isShortcut()) {
                Edge* child_1 = nullptr;
                for (int k = 0; k < m_reverse_graph[i].size(); ++k) {
                    if (m_reverse_graph[i][k].m_target == e.m_contraction_node) {
                        child_1 = &m_reverse_graph[i][k];
                        break;
                    }
                }
                Edge* child_2 = nullptr;
                for (int k = 0; k < m_reverse_graph[e.m_contraction_node].size(); ++k) {
                    if (m_reverse_graph[e.m_contraction_node][k].m_target == e.m_target) {
                        child_2 = &m_reverse_graph[e.m_contraction_node][k];
                        break;
                    }
                }
                m_reverse_graph[i][j].m_child_1 = child_1;
                m_reverse_graph[i][j].m_child_2 = child_2;
            }
        }
    }

    m_graph_contr.clear();
    m_reverse_graph_contr.clear();
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

void Graph::createCHwithoutIS(Heuristic heuristic) {
    std::cout << "Started creating CH." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty()) {
        std::cout << "Can't create CH, because graph is empty." << std::endl;
        return;
    }
    m_node_level.clear();
    m_node_level.resize(m_num_nodes);

    createContractionGraphs();

    std::vector<bool> contracted(m_num_nodes, false);
    int num_contracted = 0;

    m_contr_data.resize(1);

    m_contr_data[0] = ContractionData(m_num_nodes);
    m_contr_data[0].m_num_contracted_neighbors.clear();
    m_contr_data[0].m_num_contracted_neighbors.resize(m_num_nodes, 0);

    // initialize importance for all nodes
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        importance_pq;
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
                importance = mixedHeuristic(contracted, i, num_contracted);
                break;
            default:
                break;
        }
        importance_pq.emplace(std::make_pair(importance, i));
    }

    while (num_contracted != m_num_nodes) {
        auto contracted_node = importance_pq.top();
        importance_pq.pop();

        int new_importance = 0;
        switch (heuristic) {
            case Heuristic::IN_OUT:
                new_importance = inOutProductHeuristic(contracted, contracted_node.second);
                break;
            case Heuristic::EDGE_DIFFERENCE:
                new_importance = edgeDifferenceHeuristic(contracted, contracted_node.second);
                break;
            case Heuristic::WEIGHTED_COST:
                new_importance = weightedCostHeuristic(contracted, contracted_node.second);
                break;
            case Heuristic::MIXED:
                new_importance = mixedHeuristic(contracted, contracted_node.second, num_contracted);
                break;
            default:
                break;
        }

        while (new_importance > importance_pq.top().first) {
            importance_pq.emplace(std::make_pair(new_importance, contracted_node.second));
            contracted_node = importance_pq.top();
            importance_pq.pop();

            switch (heuristic) {
                case Heuristic::IN_OUT:
                    new_importance = inOutProductHeuristic(contracted, contracted_node.second);
                    break;
                case Heuristic::EDGE_DIFFERENCE:
                    new_importance = edgeDifferenceHeuristic(contracted, contracted_node.second);
                    break;
                case Heuristic::WEIGHTED_COST:
                    new_importance = weightedCostHeuristic(contracted, contracted_node.second);
                    break;
                case Heuristic::MIXED:
                    new_importance = mixedHeuristic(contracted, contracted_node.second, num_contracted);
                    break;
                default:
                    break;
            }
        }

        contractNode(contracted, contracted_node.second, 0);

        for (auto& fwd_shortcut : m_contr_data[0].m_shortcuts_fwd)
            m_graph_contr[fwd_shortcut.first].push_back(fwd_shortcut.second);
        for (auto& bwd_shortcut : m_contr_data[0].m_shortcuts_bwd)
            m_reverse_graph_contr[bwd_shortcut.first].push_back(bwd_shortcut.second);

        m_contr_data[0].m_shortcuts_fwd.clear();
        m_contr_data[0].m_shortcuts_bwd.clear();

        contracted[contracted_node.second] = true;
        m_node_level[contracted_node.second] = num_contracted;
        ++num_contracted;

        // std::cout << "Finished contracting: " << num_contracted << "\n";
    }

    m_contr_data.clear();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    std::cout << "Finished creating CH. Took " << elapsed.count() << " seconds" << std::endl;
}

void Graph::createCHwithIS(Heuristic heuristic) {
    std::cout << "Started creating CH." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty()) {
        std::cout << "Can't create CH, because graph is empty." << std::endl;
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
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        importance_pq;
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
        importance_pq.emplace(std::make_pair(importance, i));
    }

    while (num_contracted != m_num_nodes) {
        // create independet set
        std::vector<bool> marked(m_num_nodes, false);
        std::vector<int> independent_set;
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
            new_pq;

        // TODO: this can probably be done more efficiently
        while (!importance_pq.empty()) {
            auto contracted_node = importance_pq.top();
            importance_pq.pop();

            if (marked[contracted_node.second]) {
                new_pq.emplace(contracted_node);
                continue;
            }

            if (cur_level != 0) {
                int new_importance = 0;
                switch (heuristic) {
                    case Heuristic::IN_OUT:
                        new_importance = inOutProductHeuristic(contracted, contracted_node.second);
                        break;
                    case Heuristic::EDGE_DIFFERENCE:
                        new_importance = edgeDifferenceHeuristic(contracted, contracted_node.second);
                        break;
                    case Heuristic::WEIGHTED_COST:
                        new_importance = weightedCostHeuristic(contracted, contracted_node.second);
                        break;
                    case Heuristic::MIXED:
                        new_importance = mixedHeuristic(contracted, contracted_node.second, cur_level);
                        break;
                    default:
                        break;
                }

                while (new_importance > importance_pq.top().first) {
                    if (marked[contracted_node.second])
                        new_pq.emplace(std::make_pair(new_importance, contracted_node.second));
                    else
                        importance_pq.emplace(std::make_pair(new_importance, contracted_node.second));
                    contracted_node = importance_pq.top();
                    importance_pq.pop();

                    switch (heuristic) {
                        case Heuristic::IN_OUT:
                            new_importance = inOutProductHeuristic(contracted, contracted_node.second);
                            break;
                        case Heuristic::EDGE_DIFFERENCE:
                            new_importance = edgeDifferenceHeuristic(contracted, contracted_node.second);
                            break;
                        case Heuristic::WEIGHTED_COST:
                            new_importance = weightedCostHeuristic(contracted, contracted_node.second);
                            break;
                        case Heuristic::MIXED:
                            new_importance = mixedHeuristic(contracted, contracted_node.second, cur_level);
                            break;
                        default:
                            break;
                    }
                }
            }

            if (marked[contracted_node.second])
                new_pq.emplace(contracted_node);
            else {
                marked[contracted_node.second] = true;
                independent_set.push_back(contracted_node.second);
                for (ContractionEdge& e : m_graph_contr[contracted_node.second]) marked[e.m_target] = true;
                for (ContractionEdge& e : m_reverse_graph_contr[contracted_node.second]) marked[e.m_target] = true;
            }
        }
        importance_pq = new_pq;

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
        // std::cout << "Finished contracting: " << num_contracted << "\n";
    }

    m_contr_data.clear();
    omp_set_num_threads(m_num_threads);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    std::cout << "Finished creating CH. Took " << elapsed.count() << " seconds" << std::endl;
}

void Graph::createContractionGraphs() {
    m_graph_contr.resize(m_num_nodes);
    m_reverse_graph_contr.resize(m_num_nodes);

    for (int i = 0; i < m_num_nodes; ++i) {
        for (Edge& e : m_graph[i]) {
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
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    int num_visited_outgoing = 0;
    m_contr_data[thread_num].m_distances[start] = 0;
    m_contr_data[thread_num].m_reset_distances.push_back(start);
    pq.push(std::make_pair(0, start));

    while (!pq.empty()) {
        std::pair<int, int> cur_node = pq.top();
        pq.pop();

        if (m_contr_data[thread_num].m_distances[cur_node.second] != cur_node.first) continue;

        if (m_contr_data[thread_num].m_outgoing[cur_node.second]) ++num_visited_outgoing;
        if (cur_node.first > max_distance || num_visited_outgoing == num_outgoing) break;

        for (ContractionEdge& e : m_graph_contr[cur_node.second]) {
            if (contracted[e.m_target]) continue;
            if (m_contr_data[thread_num].m_distances[e.m_target] >
                m_contr_data[thread_num].m_distances[cur_node.second] + e.m_cost) {
                if (m_contr_data[thread_num].m_distances[e.m_target] == std::numeric_limits<int>::max())
                    m_contr_data[thread_num].m_reset_distances.push_back(e.m_target);
                m_contr_data[thread_num].m_distances[e.m_target] =
                    m_contr_data[thread_num].m_distances[cur_node.second] + e.m_cost;
                pq.push(std::make_pair(m_contr_data[thread_num].m_distances[e.m_target], e.m_target));
            }
        }
    }
}

int Graph::simplifiedHubLabelQuery(std::vector<std::pair<int, int>>& fwd_labels, int node) {
    int distance = std::numeric_limits<int>::max();
    auto fwd_iter = fwd_labels.begin();
    u_int64_t bwd_node_index = m_bwd_indices[m_node_indices[node]];
    u_int64_t bwd_next_index = m_bwd_indices[m_node_indices[node] + 1];

    while (fwd_iter != fwd_labels.end() && bwd_node_index < bwd_next_index) {
        if (fwd_iter->first == m_bwd_hub_labels[bwd_node_index].first) {
            if (fwd_iter->second + m_bwd_hub_labels[bwd_node_index].second < distance)
                distance = fwd_iter->second + m_bwd_hub_labels[bwd_node_index].second;
            ++fwd_iter;
            ++bwd_node_index;
        } else if (fwd_iter->first < m_bwd_hub_labels[bwd_node_index].first) {
            ++fwd_iter;
        } else {
            ++bwd_node_index;
        }
    }

    return distance;
}

int Graph::simplifiedHubLabelQuery(int node, std::vector<std::pair<int, int>>& bwd_labels) {
    int distance = std::numeric_limits<int>::max();
    u_int64_t fwd_node_index = m_fwd_indices[m_node_indices[node]];
    u_int64_t fwd_next_index = m_fwd_indices[m_node_indices[node] + 1];
    auto bwd_iter = bwd_labels.begin();

    while (fwd_node_index < fwd_next_index && bwd_iter != bwd_labels.end()) {
        if (bwd_iter->first == m_fwd_hub_labels[fwd_node_index].first) {
            if (bwd_iter->second + m_fwd_hub_labels[fwd_node_index].second < distance)
                distance = bwd_iter->second + m_fwd_hub_labels[fwd_node_index].second;
            ++bwd_iter;
            ++fwd_node_index;
        } else if (bwd_iter->first < m_fwd_hub_labels[fwd_node_index].first) {
            ++bwd_iter;
        } else {
            ++fwd_node_index;
        }
    }

    return distance;
}

}  // namespace labosm
