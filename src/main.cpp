#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <regex>
#include <string>
#include <vector>

#include "graph.h"
#include "graph_creator.h"
#include "server.h"

// Structure to hold benchmark results
struct BenchmarkResult {
    std::string name;
    double avg_query_time_us = 0.0;
    double avg_path_time_us = 0.0;
    double avg_pq_pops = 0.0;
    double avg_label_size = 0.0;
    double speed_up = 0.0;
    long long total_pq_pops = 0;
    long long total_query_time_us = 0;
    long long total_path_time_us = 0;
    int num_successful_queries = 0;
};

// Helper to generate CH FMI path
std::string getChFmiPath(const std::string& base_fmi_path, labosm::Heuristic heuristic) {
    std::string base_name = base_fmi_path;
    size_t pos_fmi = base_name.rfind(".fmi");
    if (pos_fmi != std::string::npos) {
        base_name.replace(pos_fmi, 4, "");
    } else {
        size_t pos_chfmi = base_name.rfind(".chfmi");
        if (pos_chfmi != std::string::npos) {
            base_name.replace(pos_chfmi, 6, "");
        }
    }

    std::string heuristic_str;
    switch (heuristic) {
        case labosm::Heuristic::IN_OUT:
            heuristic_str = "IN_OUT";
            break;
        case labosm::Heuristic::EDGE_DIFFERENCE:
            heuristic_str = "EDGE_DIFFERENCE";
            break;
        case labosm::Heuristic::WEIGHTED_COST:
            heuristic_str = "WEIGHTED_COST";
            break;
        case labosm::Heuristic::MIXED:
            heuristic_str = "MIXED";
            break;
        default:
            heuristic_str = "UNKNOWN_HEURISTIC";
    }
    return base_name + "_" + heuristic_str + ".chfmi";
}

// Function to run benchmarks
void runBenchmarks(const std::string& fmi_path, int num_queries_param, bool enable_hub_labels = true) {
    std::cout << "Starting benchmarks with " << num_queries_param << " queries on base file: " << fmi_path << std::endl;

    std::vector<std::pair<int, int>> queries(num_queries_param);
    int num_nodes_for_queries = 0;

    {
        labosm::Graph temp_graph_for_nodes(fmi_path, false, labosm::Heuristic::IN_OUT);
        num_nodes_for_queries = temp_graph_for_nodes.getNumNodes();
    }

    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<int> distrib(0, num_nodes_for_queries > 0 ? num_nodes_for_queries - 1 : 0);

    for (int i = 0; i < num_queries_param; ++i) {
        queries[i] = {distrib(rng), distrib(rng)};
    }

    std::vector<BenchmarkResult> results;
    std::vector<int> dijkstra_distances(num_queries_param);
    std::vector<bool> dijkstra_path_exists(num_queries_param);

    // Dijkstra Benchmark (Baseline)
    {
        std::cout << "Benchmarking Dijkstra (as baseline)..." << std::endl;
        BenchmarkResult dijkstra_res;
        dijkstra_res.name = "Dijkstra";

        labosm::Graph graph_dijkstra(fmi_path, false, labosm::Heuristic::IN_OUT);
        labosm::DijkstraQueryData query_data(graph_dijkstra.getNumNodes());

        for (int i = 0; i < num_queries_param; ++i) {
            query_data.reset();
            query_data.m_start = queries[i].first;
            query_data.m_end = queries[i].second;

            auto query_start_time = std::chrono::high_resolution_clock::now();
            graph_dijkstra.dijkstraQuery(query_data);
            auto query_end_time = std::chrono::high_resolution_clock::now();
            dijkstra_res.total_query_time_us +=
                std::chrono::duration_cast<std::chrono::microseconds>(query_end_time - query_start_time).count();
            dijkstra_res.total_pq_pops += query_data.num_pq_pops;

            dijkstra_distances[i] = query_data.m_distance;
            dijkstra_path_exists[i] = (query_data.m_distance != std::numeric_limits<int>::max());

            if (dijkstra_path_exists[i]) {  // Path exists
                dijkstra_res.num_successful_queries++;
                auto path_start_time = std::chrono::high_resolution_clock::now();
                graph_dijkstra.dijkstraExtractPath(query_data);
                auto path_end_time = std::chrono::high_resolution_clock::now();
                dijkstra_res.total_path_time_us +=
                    std::chrono::duration_cast<std::chrono::microseconds>(path_end_time - path_start_time).count();
            }
        }
        results.push_back(dijkstra_res);
        std::cout << "Dijkstra benchmark finished." << std::endl;
    }

    // CH and Hub Label Benchmarks (per heuristic)
    std::vector<labosm::Heuristic> ch_heuristics = {labosm::Heuristic::IN_OUT, labosm::Heuristic::EDGE_DIFFERENCE,
                                                    labosm::Heuristic::WEIGHTED_COST, labosm::Heuristic::MIXED};

    for (auto heuristic : ch_heuristics) {
        std::string heuristic_name_str;
        switch (heuristic) {
            case labosm::Heuristic::IN_OUT:
                heuristic_name_str = "IN_OUT";
                break;
            case labosm::Heuristic::EDGE_DIFFERENCE:
                heuristic_name_str = "EDGE_DIFFERENCE";
                break;
            case labosm::Heuristic::WEIGHTED_COST:
                heuristic_name_str = "WEIGHTED_COST";
                break;
            case labosm::Heuristic::MIXED:
                heuristic_name_str = "MIXED";
                break;
            default:
                heuristic_name_str = "UNKNOWN";
        }

        std::string ch_fmi_file_path = getChFmiPath(fmi_path, heuristic);

        std::cout << "Loading graph for " << heuristic_name_str << " CH/HL from " << ch_fmi_file_path << "..."
                  << std::endl;
        labosm::Graph graph_ch_hl(ch_fmi_file_path, true, 16, heuristic);
        std::cout << "  Graph for " << heuristic_name_str << " CH/HL base loaded/initialized." << std::endl;

        // --- CH Benchmark for current heuristic ---
        std::cout << "Benchmarking CH with " << heuristic_name_str << " heuristic..." << std::endl;
        BenchmarkResult ch_res;
        ch_res.name = "CH (" + heuristic_name_str + ")";
        labosm::QueryData query_data_ch(graph_ch_hl.getNumNodes());

        for (int i = 0; i < num_queries_param; ++i) {
            query_data_ch.reset();
            query_data_ch.m_start = queries[i].first;
            query_data_ch.m_end = queries[i].second;

            auto query_start_time = std::chrono::high_resolution_clock::now();
            graph_ch_hl.contractionHierarchyQuery(query_data_ch);
            auto query_end_time = std::chrono::high_resolution_clock::now();
            ch_res.total_query_time_us +=
                std::chrono::duration_cast<std::chrono::microseconds>(query_end_time - query_start_time).count();
            ch_res.total_pq_pops += query_data_ch.num_pq_pops;

            bool current_algo_path_exists = (query_data_ch.m_distance != std::numeric_limits<int>::max());

            if (current_algo_path_exists != dijkstra_path_exists[i]) {
                std::cout << "  MISMATCH (Path Existence): CH (" << heuristic_name_str << ") Query " << i << " ("
                          << queries[i].first << "->" << queries[i].second
                          << ") - CH path existence: " << current_algo_path_exists
                          << ", Dijkstra path existence: " << dijkstra_path_exists[i] << std::endl;
            } else if (current_algo_path_exists && query_data_ch.m_distance != dijkstra_distances[i]) {
                std::cout << "  MISMATCH (Distance): CH (" << heuristic_name_str << ") Query " << i << " ("
                          << queries[i].first << "->" << queries[i].second
                          << ") - CH dist: " << query_data_ch.m_distance << ", Dijkstra dist: " << dijkstra_distances[i]
                          << std::endl;
            }

            if (current_algo_path_exists) {
                ch_res.num_successful_queries++;
                auto path_start_time = std::chrono::high_resolution_clock::now();
                graph_ch_hl.contractionHierarchyExtractPath(query_data_ch);
                auto path_end_time = std::chrono::high_resolution_clock::now();
                ch_res.total_path_time_us +=
                    std::chrono::duration_cast<std::chrono::microseconds>(path_end_time - path_start_time).count();
            }
        }
        results.push_back(ch_res);
        std::cout << "CH with " << heuristic_name_str << " benchmark finished." << std::endl;

        // --- HL Benchmark for current heuristic (using the same graph_ch_hl) ---
        if (enable_hub_labels) {
            std::cout << "Benchmarking Hub Labels on " << heuristic_name_str << " CH base..." << std::endl;
            BenchmarkResult hl_res;
            hl_res.name = "Hub Labels (" + heuristic_name_str + " CH)";
            std::cout << "  Creating Hub Labels for " << heuristic_name_str << " base..." << std::endl;
            auto hl_create_start = std::chrono::high_resolution_clock::now();
            graph_ch_hl.createHubLabelsWithIS();
            auto hl_create_end = std::chrono::high_resolution_clock::now();
            long long create_time_ms =
                std::chrono::duration_cast<std::chrono::milliseconds>(hl_create_end - hl_create_start).count();
            std::cout << "  Hub Labels for " << heuristic_name_str << " created in " << create_time_ms << " ms."
                      << std::endl;
            labosm::QueryData query_data_hl(graph_ch_hl.getNumNodes());
            for (int i = 0; i < num_queries_param; ++i) {
                query_data_hl.reset();
                query_data_hl.m_start = queries[i].first;
                query_data_hl.m_end = queries[i].second;
                std::pair<uint32_t, uint32_t> hub_indices;
                auto query_start_time = std::chrono::high_resolution_clock::now();
                hub_indices = graph_ch_hl.hubLabelQuery(query_data_hl);
                auto query_end_time = std::chrono::high_resolution_clock::now();
                hl_res.total_query_time_us +=
                    std::chrono::duration_cast<std::chrono::microseconds>(query_end_time - query_start_time).count();
                hl_res.total_pq_pops += query_data_hl.num_pq_pops;
                bool current_algo_path_exists_hl = (query_data_hl.m_distance != std::numeric_limits<int>::max());
                if (current_algo_path_exists_hl != dijkstra_path_exists[i]) {
                    std::cout << "  MISMATCH (Path Existence): HL (" << heuristic_name_str << " CH) Query " << i << " ("
                              << queries[i].first << "->" << queries[i].second
                              << ") - HL path existence: " << current_algo_path_exists_hl
                              << ", Dijkstra path existence: " << dijkstra_path_exists[i] << std::endl;
                } else if (current_algo_path_exists_hl && query_data_hl.m_distance != dijkstra_distances[i]) {
                    std::cout << "  MISMATCH (Distance): HL (" << heuristic_name_str << " CH) Query " << i << " ("
                              << queries[i].first << "->" << queries[i].second
                              << ") - HL dist: " << query_data_hl.m_distance
                              << ", Dijkstra dist: " << dijkstra_distances[i] << std::endl;
                }
                if (current_algo_path_exists_hl) {
                    hl_res.num_successful_queries++;
                    auto path_start_time = std::chrono::high_resolution_clock::now();
                    graph_ch_hl.hubLabelExtractPath(query_data_hl, hub_indices);
                    auto path_end_time = std::chrono::high_resolution_clock::now();
                    hl_res.total_path_time_us +=
                        std::chrono::duration_cast<std::chrono::microseconds>(path_end_time - path_start_time).count();
                }
            }
            hl_res.avg_label_size = graph_ch_hl.averageLabelSize();
            results.push_back(hl_res);
            std::cout << "  Clearing Hub Labels for " << heuristic_name_str << " base..." << std::endl;
            graph_ch_hl.clearHubLabel();
            std::cout << "Hub Labels on " << heuristic_name_str << " CH benchmark finished." << std::endl;
        } else {
            std::cout << "Skipping Hub Label benchmarks for " << heuristic_name_str << " due to flag.\n";
        }
    }

    // Calculate averages and print results
    std::cout << "\n--- Benchmark Results (" << num_queries_param << " queries each) ---\n" << std::endl;
    std::cout << "| Algorithm                 | Avg. Query Time (us) | Avg. Path Time (us) | Avg. PQ Pops | Avg. Label "
                 "Size | Speed-up      |\n";
    std::cout << "|---------------------------|----------------------|---------------------|--------------|------------"
                 "-----|---------------|\n";

    double dijkstra_avg_query_time = 0.0;
    if (!results.empty() && results[0].name == "Dijkstra" && num_queries_param > 0) {
        dijkstra_avg_query_time = static_cast<double>(results[0].total_query_time_us) / num_queries_param;
    }

    for (auto& res : results) {
        if (num_queries_param > 0) {
            res.avg_query_time_us = static_cast<double>(res.total_query_time_us) / num_queries_param;
            res.avg_pq_pops = static_cast<double>(res.total_pq_pops) / num_queries_param;
        }
        if (res.num_successful_queries > 0) {
            res.avg_path_time_us = static_cast<double>(res.total_path_time_us) / res.num_successful_queries;
        } else {
            res.avg_path_time_us = 0;
        }

        if (dijkstra_avg_query_time > 0 && res.avg_query_time_us > 0) {
            res.speed_up = dijkstra_avg_query_time / res.avg_query_time_us;
        } else {
            res.speed_up = 0.0;
        }

        std::cout << "| " << res.name << " | " << res.avg_query_time_us << " | " << res.avg_path_time_us << " | "
                  << res.avg_pq_pops;
        if (res.name.find("Hub Labels") != std::string::npos) {
            std::cout << " | " << res.avg_label_size;
        } else {
            std::cout << " | N/A";
        }
        std::cout << " | " << res.speed_up << "x";
        std::cout << " |\n";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc <= 2) {
        std::cout << "Usage ./labosm <mode> <args>\n" << "\n";
        std::cout
            << "Supported modes: generate_points points_to_fmi simpleserver advancedserver create_chfmi benchmark\n"
            << "\n";

        std::cout << "Example: ./labosm generate_points input.osm.pbf output_prefix num_points" << "\n";
        std::cout << "Example: ./labosm generate_points input.osm.pbf output_prefix num_points image_path\n";
        std::cout << "Generates num_points random points on the sphere, filters them based on the given coastlines and "
                     "writes them to a geojson file\n";
        std::cout << "Creates the files output_prefix_coastlines.geojson and output_prefix_filtered_points.geojson\n";
        std::cout << "If image_path is given, the points are filtered based on the image\n" << "\n";

        std::cout << "Example: ./labosm points_to_fmi filtered_points.geojson output.fmi\n" << "\n";
        std::cout << "Converts the filtered points to a fmi file\n" << "\n";

        std::cout << "Example: ./labosm create_chfmi input.fmi heuristic_name\n";
        std::cout << "Takes an FMI file, generates Contraction Hierarchy data using the specified heuristic, and saves "
                     "it to input_heuristicName.chfmi\n";
        std::cout << "Supported heuristics: IN_OUT, EDGE_DIFFERENCE, WEIGHTED_COST, MIXED\n" << "\n";

        std::cout << "Example: ./labosm simpleserver input.fmi" << "\n";
        std::cout << "Simple Server only supports dijkstra for routing" << "\n";
        std::cout << "Example: ./labosm advancedserver input.fmi" << "\n";
        std::cerr << "Advanced Server supports CH and Hub Labeling for routing" << "\n";

        std::cout << "Example: ./labosm benchmark input.fmi [num_queries]" << "\n";
        std::cout << "Runs benchmarks with num_queries (default 10000) on the given fmi file and shows the results"
                  << "\n";

        return 1;
    }

    if (argv[1] == std::string("simpleserver")) {
        if (argc != 3) {
            std::cerr << "Usage: ./labosm simpleserver input.fmi" << "\n";
            return 1;
        }
        labosm::server::simpleServer(argv[2]);
        return 0;
    } else if (argv[1] == std::string("advancedserver")) {
        if (argc < 3) {
            std::cerr << "Usage: ./labosm advancedserver input.fmi [--no-hub-labels]" << "\n";
            return 1;
        }
        std::string fmi_path = argv[2];
        bool enable_hub_labels = true;
        for (int i = 3; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--no-hub-labels") {
                enable_hub_labels = false;
            }
        }
        labosm::server::advancedServer(fmi_path, enable_hub_labels);
        return 0;
    } else if (argv[1] == std::string("generate_points")) {
        if (argc < 5 || argc > 6) {
            std::cerr << "Usage: ./labosm generate_points input.osm.pbf output_prefix num_points [image_path]" << "\n";
            return 1;
        }
        labosm::GraphCreator graph_creator;
        std::string image_path = (argc == 6) ? argv[5] : "";
        graph_creator.generatePointsAndFilter(argv[2], std::stoi(argv[4]), argv[3], image_path != "", image_path);
        return 0;
    } else if (argv[1] == std::string("points_to_fmi")) {
        if (argc != 4) {
            std::cerr << "Usage: ./labosm points_to_fmi filtered_points.geojson output.fmi" << "\n";
            return 1;
        }
        labosm::GraphCreator graph_creator;
        graph_creator.generateGraph(argv[2], argv[3]);
        return 0;
    } else if (argv[1] == std::string("create_chfmi")) {
        if (argc != 4) {  // Expect 4 arguments: ./labosm create_chfmi input.fmi heuristic_name
            std::cerr << "Usage: ./labosm create_chfmi input.fmi heuristic_name" << "\n";
            std::cerr << "Supported heuristics: IN_OUT, EDGE_DIFFERENCE, WEIGHTED_COST, MIXED" << "\n";
            return 1;
        }

        std::string heuristic_arg = argv[3];
        labosm::Heuristic ch_heuristic_to_use;

        if (heuristic_arg == "IN_OUT") {
            ch_heuristic_to_use = labosm::Heuristic::IN_OUT;
        } else if (heuristic_arg == "EDGE_DIFFERENCE") {
            ch_heuristic_to_use = labosm::Heuristic::EDGE_DIFFERENCE;
        } else if (heuristic_arg == "WEIGHTED_COST") {
            ch_heuristic_to_use = labosm::Heuristic::WEIGHTED_COST;
        } else if (heuristic_arg == "MIXED") {
            ch_heuristic_to_use = labosm::Heuristic::MIXED;
        } else {
            std::cerr << "Unknown heuristic: " << heuristic_arg << "\n";
            std::cerr << "Supported heuristics: IN_OUT, EDGE_DIFFERENCE, WEIGHTED_COST, MIXED" << "\n";
            return 1;
        }

        // Instantiate Graph object with CH enabled and the chosen heuristic.
        // The Graph constructor will handle CH creation and writing the .chfmi file.
        // Using the constructor that supports IS (num_threads = 16 as per previous user preference).
        labosm::Graph graph(argv[2], true, 16, ch_heuristic_to_use);
        std::cout << "Successfully created .chfmi file." << std::endl;
        return 0;
    } else if (argv[1] == std::string("benchmark")) {
        if (argc < 3) {  // Expect at least 3 arguments
            std::cerr << "Usage: ./labosm benchmark input.fmi [num_queries] [--no-hub-labels]" << "\n";
            return 1;
        }
        std::string fmi_path = argv[2];
        int num_queries = 10000;
        bool enable_hub_labels = true;
        for (int i = 3; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--no-hub-labels") {
                enable_hub_labels = false;
            } else {
                try {
                    num_queries = std::stoi(arg);
                } catch (...) {
                    std::cerr << "Invalid argument: " << arg << "\n";
                    return 1;
                }
            }
        }
        runBenchmarks(fmi_path, num_queries, enable_hub_labels);
        return 0;
    } else {
        std::cerr << "Unknown mode: " << argv[1] << "\n";
        return 1;
    }

    return 0;
}
