#pragma once

#include <omp.h>
#include <stdint.h>

#include <cstdint>
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
     * @param data that is reused for multiple queries
     */
    void dijkstraQuery(DijkstraQueryData& data);

    /**
     * @brief Extracts the path from the data structure.
     *  Make sure to call dijkstraQuery before using this function.
     *
     * @param data to extract path from
     */
    void dijkstraExtractPath(DijkstraQueryData& data);

    /**
     * @brief Works similarly to bidirectionalDijkstraQuery.
     *  Additionally uses CH information to reduce query time.
     *
     * @param data set m_path_needed = true to make path extractable
     */
    void contractionHierarchyQuery(QueryData& data);

    /**
     * @brief Extracts the path from the data structure.
     *  Make sure to call contractionHierarchyQuery before using this function.
     *
     * @param data to extract path from
     */
    void contractionHierarchyExtractPath(QueryData& data);

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
    std::pair<uint64_t, uint64_t> hubLabelQuery(QueryData& data);

    /**
     * @brief Extracts the path from the data structure.
     *  Make sure to call hubLabelQuery before using this function.
     *
     * @param data to extract path from
     * @param hub_indices indices of the hub labels
     */
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
        std::vector<uint64_t>().swap(m_fwd_indices);
        std::vector<uint64_t>().swap(m_bwd_indices);
        std::vector<std::tuple<uint32_t, uint32_t, uint16_t, uint16_t>>().swap(m_fwd_hub_labels);
        std::vector<std::tuple<uint32_t, uint32_t, uint16_t, uint16_t>>().swap(m_bwd_hub_labels);
    }

    /**
     * @brief Clears all the dijkstra graph data from memory.
     * Useful for freeing memory when dijkstra functionality is no longer needed.
     */
    void clearDijkstraGraph() {
        std::vector<uint64_t>().swap(m_dijkstra_indices);
        std::vector<SimpleEdge>().swap(m_dijkstra_edges);
    }

   private:
    bool m_ch_available;  // ch/CH = Contraction Hierarchy
    int m_num_nodes;
    bool m_is;          // determines whether IS are used
    int m_num_threads;  // sets threads when using independent sets (IS)

    // dijkstra graph data structure with flat edge storage
    std::vector<uint64_t> m_dijkstra_indices;  // indices for each node to the corresponding edges
    std::vector<SimpleEdge> m_dijkstra_edges;  // flat vector containing all edges

    // used for CH/HL
    std::vector<std::vector<Edge>> m_graph;
    std::vector<std::vector<Edge>> m_reverse_graph;
    std::vector<int> m_node_level;
    std::vector<std::pair<double, double>> m_node_coords;

    // NOTE: This memory is freed once the contraction hierarchy is created
    // only used to share state when calculating CH
    std::vector<ContractionData> m_contr_data;
    // copies of the graph with additional data, only used when creating CH
    std::vector<std::vector<ContractionEdge>> m_graph_contr;
    std::vector<std::vector<ContractionEdge>> m_reverse_graph_contr;

    // needed for hub label creation since, we go from the highest level to the lowest level
    std::vector<int> m_level_indices_sorted;
    std::vector<int> m_node_indices;
    // indices for each nodes to the corresponding hub labels
    std::vector<uint64_t> m_fwd_indices;
    std::vector<uint64_t> m_bwd_indices;
    // first is hub node, second is distance to the hub node, third is edge_index in the graph (edge to that hub
    // node), fourth is the offset to the original hub within the node's hub labels
    // Using smaller types: uint32_t for node IDs and distances, uint16_t for edge indices and hub label offsets
    std::vector<std::tuple<uint32_t, uint32_t, uint16_t, uint16_t>> m_fwd_hub_labels;
    std::vector<std::tuple<uint32_t, uint32_t, uint16_t, uint16_t>> m_bwd_hub_labels;

    /**
     * @brief Reads the graph from a file.
     * Should be a fmi file.
     * @param path
     */
    void readGraph(const std::string& path);

    /**
     * @brief Reads a simple graph from a file for dijkstra queries.
     * @param path
     */
    void readSimpleGraph(const std::string& path);

    /**
     * @brief Writes the graph to a file in the FMI format including CH information.
     *
     * @param filename The name of the file to write to.
     */
    void writeGraphToCHFMI(const std::string& filename);

    /**
     * @brief Reads the graph from a CHFMI file.
     * Should be a chfmi file.
     * @param path
     */
    void readGraphFromCHFMI(const std::string& path);

    /**
     * @brief Creates the reverse graph.
     * Depending on the mode it calls the normal or CH version.
     */
    void createReverseGraph();

    /**
     * @brief If CH is used, this function creates the reverse graph.
     * The forward graph only stores upward edges and the reverse graph only stores downward edges.
     * This saves memory and we don't need to make an upward/downward edge test during the query.
     */
    void createReverseGraphCH();

    /**
     * @brief Creates the reverse graph normally.
     * So edge edges are stored in both directions.
     * Is used when using the simpleserver/dijkstra mode.
     */
    void createReverseGraphNormal();

    /**
     * @brief Calculates hierarchy and adds shortcuts to graph without using independent sets.
     * This means that many nodes have the same level/hierarchy value.
     *
     * @param heuristic
     */
    void createCHwithIS(Heuristic heuristic);

    /**
     * @brief Creates intermediate contraction graphs.
     * Are saved in m_graph_contr and m_reverse_graph_contr.
     * Are just used for the contraction hierarchy creation.
     */
    void createContractionGraphs();

    /**
     * @brief Just counts the number of non-contracted incoming and outgoing edges and multiplies them.
     * This is the baseline heuristic.
     * It is very fast, but not recommended for larger graphs.
     */
    int inOutProductHeuristic(std::vector<bool>& contracted, int node);

    /**
     * @brief The classical edge difference heuristic.
     * Counts the number deleted edges and the number of added shortcuts and takes the difference.
     */
    int edgeDifferenceHeuristic(std::vector<bool>& contracted, int node);

    /**
     * @brief Heuristic that I developed during my bachelor thesis.
     * Finds the largest newly added shortcut and weights it together with the inout product heuristic.
     * It performs the best when not using IS.
     */
    int weightedCostHeuristic(std::vector<bool>& contracted, int node);

    /**
     * @brief It mixes my weighted cost heuristic with the heuristic found in the paper A Hub-Based Labeling Algorithm
     * for Shortest Paths on Road Networks by I. Abraham et al.
     * (https://www.microsoft.com/en-us/research/wp-content/uploads/2010/12/HL-TR.pdf).
     * It performs the best when using IS.
     *
     * It takes the largest newly added shortcut, edge difference, number of contracted neighbors, number of underlying
     * edges and the level that would be assigned to the node into account.
     */
    int mixedHeuristic(std::vector<bool>& contracted, int node, int cur_level);

    /**
     * @brief Performs the actual contraction of the node.
     * It adds the shortcuts to the graph and "removes/deactivates" the node from the graph.
     */
    void contractNode(std::vector<bool>& contracted, int contracted_node, int thread_num);

    /**
     * @brief A simpliefed dijkstra that stops much earlier to check if need to create a shortcut between the neighbours
     * of the currentlcy contracted node.
     */
    void contractionDijkstra(int start, int contracted_node, std::vector<bool>& contracted, int num_outgoing,
                             int max_distance, int thread_num);

    /**
     * @brief Simplified hub label query that is used during the construction of the hub labels.
     * Is used for the fwd hub label creation.
     * @param fwd_labels
     * @param node
     * @return distance to the node
     */
    int simplifiedHubLabelQuery(std::vector<std::tuple<uint32_t, uint32_t, uint16_t, uint16_t>>& fwd_labels, int node);

    /**
     * @brief Simplified hub label query that is used during the construction of the hub labels.
     * Is used for the bwd hub label creation.
     * @param node
     * @param bwd_labels
     * @return distance to the node
     */
    int simplifiedHubLabelQuery(int node, std::vector<std::tuple<uint32_t, uint32_t, uint16_t, uint16_t>>& bwd_labels);

    /**
     * @brief Unpacks shortcut edges into the original edges.
     * Is used to extract the original path from contraction query and hub label query.
     */
    void unpackEdge(const std::tuple<int, int, bool>& edge, std::vector<int>& path);
};

}  // namespace labosm
