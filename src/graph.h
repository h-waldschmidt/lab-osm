#pragma once

#include <omp.h>
#include <stdint.h>

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "helper.h"
namespace labosm {

/**
 * @brief Stores graph data structure and allows for various operations on graph data structure.
 *
 */
class Graph {
   public:
    /**
     * @brief Construct a new Graph object independent sets (IS)
     *
     * @param path          specifies location of graph
     * @param read_mode     specifies if CH is pre-calculated and stored in graph-file
     * @param ch_available  specifies if CH should be used
     * @param ch_heuristic  specifies which heuristic to use
     * @param dist_mode     specifies whether to use travel time or distance in meter metric
     */
    Graph(const std::string& path, bool ch_available, Heuristic ch_heuristic);

    // constructor with independent sets (IS)
    /**
     * @brief Construct a new Graph object with independent sets (IS)
     *
     * @param path          specifies location of graph
     * @param read_mode     specifies if CH is pre-calculated and stored in graph-file
     * @param ch_available  specifies if CH should be used
     * @param num_threads   specifies the number of threads to use in createHubLabelsWithIS and createCHwithIS functions
     * @param ch_heuristic  specifies which heuristic to use
     * @param dist_mode     specifies whether to use travel time or distance in meter metric
     */
    Graph(const std::string& path, bool ch_available, int num_threads, Heuristic ch_heuristic);

    Graph() = default;
    ~Graph() = default;

    /**
     * @brief Conventional implementation of the dijkstra algorithm
     * Should only be used to compare query times to other methods.
     *
     * @param graph since the function is static the graph vector must be passed
     * @param start node
     * @param end node
     * @return distances between start and end
     */
    void dijkstraQuery(DijkstraQueryData& data);

    void dijkstraExtractPath(DijkstraQueryData& data);

    /**
     * @brief Extracts path from QueryData.
     * dijkstraQuery function with m_path_needed = true needs to be called beforehand.
     *
     * @param data
     */
    void bidirectionalDijkstraGetPath(QueryData& data);

    /**
     * @brief Works similarly to bidirectionalDijkstraQuery.
     *  Additionally uses CH information to reduce query time.
     *
     * @param data set m_path_needed = true to make path extractable
     */
    void contractionHierarchyQuery(QueryData& data);

    void contractionHierarchyExtractPath(QueryData& data);

    /**
     * @brief Create a hub labeling when not using IS for CH.
     */
    void createHubLabelsWithoutIS();

    /**
     * @brief Create a hub labeling when using IS for CH.
     * Implementation uses OpenMP for parallelization with specified number of threads.
     */
    void createHubLabelsWithIS();

    /**
     * @brief Calculates distance using the created hub labeling.
     * Make sure hub labeling has been created before using this function.
     *
     * @param data
     */
    std::pair<uint32_t, uint32_t> hubLabelQuery(QueryData& data);

    void hubLabelExtractPath(QueryData& data, std::pair<int, int> hub_indices);

    /**
     * @brief Calculates the distance of two given points based on the Haversine Formula
     *
     * @param lat_1
     * @param lon_1
     * @param lat_2
     * @param lon_2
     * @return distance between the points in meter
     */
    static int greatCircleDistance(double lat_1, double lon_1, double lat_2, double lon_2);

    double averageLabelSize();

    int maxLabelSize();

    int getNearestNode(double latitude, double longitude);

    int getNumNodes() { return m_num_nodes; }

    std::pair<double, double> getNodeCoords(int node) { return m_node_coords[node]; }

    void setNumThreads(int num_of_threads) { m_num_threads = num_of_threads; }

    /**
     * @brief Clears all the hub label data from memory.
     * Useful when ESC has been extracted from hub labeling and hub labeling is no longer required.
     *
     */
    void clearHubLabel() {
        std::vector<uint32_t>().swap(m_fwd_indices);
        std::vector<uint32_t>().swap(m_bwd_indices);
        std::vector<std::tuple<int, int, int, int>>().swap(m_fwd_hub_labels);
        std::vector<std::tuple<int, int, int, int>>().swap(m_bwd_hub_labels);
    }

   private:
    bool m_ch_available;  // ch/CH = Contraction Hierarchy
    int m_num_nodes;
    bool m_is;          // determines whether IS are used
    int m_num_threads;  // sets threads when using independent sets (IS)

    std::vector<std::vector<Edge>> m_graph;
    std::vector<std::vector<Edge>> m_reverse_graph;
    std::vector<int> m_node_level;
    std::vector<std::pair<double, double>> m_node_coords;

    // only used to share state when calculating CH
    std::vector<ContractionData> m_contr_data;
    // copies of the graph with additional data, only used when creating CH
    std::vector<std::vector<ContractionEdge>> m_graph_contr;
    std::vector<std::vector<ContractionEdge>> m_reverse_graph_contr;

    std::vector<int> m_level_indices_sorted;
    std::vector<int> m_node_indices;
    std::vector<uint32_t> m_fwd_indices;
    std::vector<uint32_t> m_bwd_indices;
    // first is hub node, second is distance to the hub node, third is edge_index in the graph (edge to that hub
    // node), fourth is the index to the original hub
    std::vector<std::tuple<int, int, int, int>> m_fwd_hub_labels;
    std::vector<std::tuple<int, int, int, int>> m_bwd_hub_labels;

    void readGraph(const std::string& path);

    void createReverseGraph();

    void createReverseGraphCH();

    void createReverseGraphNormal();

    /**
     * @brief Calculates hierarchy and adds shortcuts to graph without using independent sets.
     * This means that each node has a distinct level/hierarchy value.
     *
     * @param heuristic
     */
    void createCHwithoutIS(Heuristic heuristic);

    /**
     * @brief Calculates hierarchy and adds shortcuts to graph without using independent sets.
     * This means that many nodes have the same level/hierarchy value.
     *
     * @param heuristic
     */
    void createCHwithIS(Heuristic heuristic);

    void createContractionGraphs();

    int inOutProductHeuristic(std::vector<bool>& contracted, int node);

    int edgeDifferenceHeuristic(std::vector<bool>& contracted, int node);

    int weightedCostHeuristic(std::vector<bool>& contracted, int node);

    int mixedHeuristic(std::vector<bool>& contracted, int node, int cur_level);

    void contractNode(std::vector<bool>& contracted, int contracted_node, int thread_num);

    void contractionDijkstra(int start, int contracted_node, std::vector<bool>& contracted, int num_outgoing,
                             int max_distance, int thread_num);

    int simplifiedHubLabelQuery(std::vector<std::tuple<int, int, int, int>>& fwd_labels, int node);

    int simplifiedHubLabelQuery(int node, std::vector<std::tuple<int, int, int, int>>& bwd_labels);

    void unpackEdge(const std::tuple<int, int, bool>& edge, std::vector<int>& path);
};

}  // namespace labosm
