#include <cstdint>
#include <fstream>
#include <iostream>
#include <osmium/handler.hpp>
#include <osmium/handler/node_locations_for_ways.hpp>
#include <osmium/index/map/sparse_mem_array.hpp>
#include <osmium/io/any_input.hpp>
#include <osmium/osm/location.hpp>
#include <osmium/visitor.hpp>
#include <regex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "graph.h"
#include "helper.h"

using NodeMap = std::unordered_map<uint64_t, std::pair<double, double>>;
using WayList = std::unordered_map<uint64_t, std::vector<uint64_t>>;

NodeMap coastline_nodes;
WayList coastline_ways;

int parseLine(std::string line) {
    return stoi(std::regex_replace(line, std::regex("[^0-9]*([0-9]+).*"), std::string("$1")));
}

// https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
int getMemoryUsage() {  // Note: this value is in KB!
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

class CoastlineHandler : public osmium::handler::Handler {
   public:
    void way(const osmium::Way& way) {
        if (way.tags().has_tag("natural", "coastline")) {
            std::vector<uint64_t> way_nodes;
            for (const auto& node_ref : way.nodes()) {
                osmium::Location loc = node_ref.location();
                if (loc.valid()) {
                    uint64_t id = node_ref.ref();
                    way_nodes.push_back(id);
                    coastline_nodes[id] = {loc.lon(), loc.lat()};
                }
            }
            if (!way_nodes.empty()) {
                coastline_ways[way_nodes.front()] = way_nodes;
            }
        }
    }
};

void write_geojson(const std::string& output_file, const NodeMap& nodes, const WayList& ways) {
    std::ofstream out(output_file);
    out << R"({"type": "FeatureCollection", "features": [)";

    bool first = true;
    for (const auto& way : ways) {
        if (!first) out << ",";
        first = false;

        out << R"({"type": "Feature","geometry":{"type": "LineString","coordinates":[)";
        const auto& way_nodes = way.second;
        for (size_t i = 0; i < way_nodes.size(); ++i) {
            uint64_t id = way_nodes[i];
            auto it = nodes.find(id);
            if (it != nodes.end()) {
                const auto& [lon, lat] = it->second;
                if (i > 0) out << ",";
                out << "[" << lon << "," << lat << "]";
            }
        }
        out << R"(]},"properties":{}})";
    }

    out << "]}" << '\n';
}

void merge_ways() {
    size_t old_size = coastline_ways.size();
    while (true) {
        for (auto it1 = coastline_ways.begin(); it1 != coastline_ways.end(); ++it1) {
            uint64_t last_id = it1->second.back();
            if (last_id == it1->first) continue;

            auto it2 = coastline_ways.find(last_id);
            if (it2 != coastline_ways.end()) {
                it1->second.insert(it1->second.end(), it2->second.begin(), it2->second.end());
                coastline_ways.erase(it2);
            }
        }

        if (coastline_ways.size() == old_size || coastline_ways.size() == 1) {
            break;
        }

        old_size = coastline_ways.size();
    }
}
/*
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: extract_coastlines input.osm.pbf output.geojson\n";
        return 1;
    }

    const char* input_filename = argv[1];
    const char* output_filename = argv[2];

    auto start_reading = std::chrono::steady_clock::now();

    osmium::io::Reader reader{input_filename, osmium::osm_entity_bits::node | osmium::osm_entity_bits::way};

    using index_type = osmium::index::map::SparseMemArray<osmium::unsigned_object_id_type, osmium::Location>;
    index_type index;
    osmium::handler::NodeLocationsForWays<index_type> location_handler{index};
    location_handler.ignore_errors();

    CoastlineHandler handler;

    osmium::apply(reader, location_handler, handler);
    reader.close();
    auto end_reading = std::chrono::steady_clock::now();

    std::cout << "Reading time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_reading - start_reading).count() << " ms\n";

    auto start_merging = std::chrono::steady_clock::now();
    merge_ways();
    auto end_merging = std::chrono::steady_clock::now();
    std::cout << "Merging time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_merging - start_merging).count() << " ms\n";

    // find the longest array
    size_t max_size = 0;
    std::vector<uint64_t> longest_way;
    for (const auto& way : coastline_ways) {
        if (way.second.size() > max_size) {
            max_size = way.second.size();
            longest_way = way.second;
        }
    }
    coastline_ways.clear();
    coastline_ways[longest_way.front()] = longest_way;

    write_geojson(output_filename, coastline_nodes, coastline_ways);

    return 0;
}
*/

#include "httplib.h"

void simpleServer(const std::string& fmi_file) {
    labosm::Graph g(fmi_file, false, 1, labosm::Heuristic::IN_OUT);

    labosm::DijkstraQueryData data(g.getNumNodes());

    httplib::Server svr;
    // first the static content: website etc.
    svr.set_mount_point("/", "static");

    svr.Get(R"(/api/)", [&](const httplib::Request& req, httplib::Response& res) {
        res.set_content(R"({"type": "simple"})", "application/json");
    });

    // get the nearest node for given coords
    svr.Get(R"(/api/nearest_node)", [&](const httplib::Request& req, httplib::Response& res) {
        if (req.has_param("lat") && req.has_param("lon")) {
            double lat = std::stod(req.get_param_value("lat"));
            double lon = std::stod(req.get_param_value("lon"));
            int nearest_node = g.getNearestNode(lat, lon);
            std::pair<double, double> coords = g.getNodeCoords(nearest_node);
            res.status = 200;
            res.set_content(R"({"node": )" + std::to_string(nearest_node) + R"(,"lat": )" +
                                std::to_string(coords.first) + R"(,"lon": )" + std::to_string(coords.second) + R"(})",
                            "application/json");
        } else {
            res.status = 400;
            res.set_content(R"({"error": "Missing lat or lon parameter"})", "application/json");
        }
    });

    // dijkstra query
    svr.Get(R"(/api/dijkstra)", [&](const httplib::Request& req, httplib::Response& res) {
        if (req.has_param("start") && req.has_param("end")) {
            int start = std::stoi(req.get_param_value("start"));
            int end = std::stoi(req.get_param_value("end"));
            data.m_start = start;
            data.m_end = end;
            g.dijkstraQuery(data);
            g.dijkstraExtractPath(data);
            res.status = 200;

            std::vector<std::pair<double, double>> coords;
            for (int node : data.m_path) {
                coords.push_back(g.getNodeCoords(node));
            }
            std::string path_str = "{\"type\": \"LineString\", \"coordinates\": [";
            for (size_t i = 0; i < coords.size(); ++i) {
                if (i > 0) path_str += ",";
                path_str += "[" + std::to_string(coords[i].first) + "," + std::to_string(coords[i].second) + "]";
            }
            path_str += "]}";

            std::cout << "Distance: " << data.m_distance << std::endl;

            res.set_content(R"({"distance": )" + std::to_string(data.m_distance) + R"(,"path": )" + path_str + R"(})",
                            "application/json");
        } else {
            res.status = 400;
            res.set_content(R"({"error": "Missing start or end parameter"})", "application/json");
        }
    });

    std::cout << "Server started on port 8080" << std::endl;
    svr.listen("0.0.0.0", 8080);
}

void advancedServer(const std::string& fmi_file) {
    labosm::Graph g(fmi_file, true, 4, labosm::Heuristic::IN_OUT);
    labosm::QueryData data(g.getNumNodes());
    httplib::Server svr;
    // first the static content: website etc.
    svr.set_mount_point("/", "static");

    svr.Get(R"(/api/)", [&](const httplib::Request& req, httplib::Response& res) {
        res.set_content(R"({"type": "complex"})", "application/json");
    });

    // get the nearest node for given coords
    svr.Get(R"(/api/nearest_node)", [&](const httplib::Request& req, httplib::Response& res) {
        if (req.has_param("lat") && req.has_param("lon")) {
            double lat = std::stod(req.get_param_value("lat"));
            double lon = std::stod(req.get_param_value("lon"));
            int nearest_node = g.getNearestNode(lat, lon);
            std::pair<double, double> coords = g.getNodeCoords(nearest_node);
            res.status = 200;
            res.set_content(R"({"node": )" + std::to_string(nearest_node) + R"(,"lat": )" +
                                std::to_string(coords.first) + R"(,"lon": )" + std::to_string(coords.second) + R"(})",
                            "application/json");
        } else {
            res.status = 400;
            res.set_content(R"({"error": "Missing lat or lon parameter"})", "application/json");
        }
    });

    // CH query
    svr.Get(R"(/api/ch)", [&](const httplib::Request& req, httplib::Response& res) {
        if (req.has_param("start") && req.has_param("end")) {
            int start = std::stoi(req.get_param_value("start"));
            int end = std::stoi(req.get_param_value("end"));
            data.m_start = start;
            data.m_end = end;
            auto start_time = std::chrono::steady_clock::now();
            g.contractionHierarchyQuery(data);
            auto end_time = std::chrono::steady_clock::now();
            std::cout << "CH Query Time: "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()
                      << std::endl;

            start_time = std::chrono::steady_clock::now();
            g.contractionHierarchyExtractPath(data);
            end_time = std::chrono::steady_clock::now();
            std::cout << "CH Path Extraction Time: "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()
                      << std::endl;

            res.status = 200;

            std::vector<std::pair<double, double>> coords;
            for (int node : data.m_shortest_path) {
                coords.push_back(g.getNodeCoords(node));
            }

            std::cout << "Start: " << data.m_start << std::endl;
            std::cout << "End: " << data.m_end << std::endl;
            std::cout << "Distance: " << data.m_distance << std::endl;
            std::cout << "Meeting Node: " << data.m_meeting_node << std::endl;

            std::string path_str = "{\"type\": \"LineString\", \"coordinates\": [";
            for (size_t i = 0; i < coords.size(); ++i) {
                if (i > 0) path_str += ",";
                path_str += "[" + std::to_string(coords[i].first) + "," + std::to_string(coords[i].second) + "]";
            }
            path_str += "]}";

            std::cout << "Distance: " << data.m_distance << std::endl;

            res.set_content(R"({"distance": )" + std::to_string(data.m_distance) + R"(,"path": )" + path_str + R"(})",
                            "application/json");
        } else {
            res.status = 400;
            res.set_content(R"({"error": "Missing start or end parameter"})", "application/json");
        }
    });

    std::cout << "Server started on port 8080" << std::endl;
    svr.listen("0.0.0.0", 8080);
}

int main(int argc, char* argv[]) {
    if (argc <= 2) {
        std::cout << "Usage ./labosm <mode> <args>\n" << std::endl;
        std::cout << "Supported modes: osmtogeojson osmtofmi simpleserver advancedserver\n" << std::endl;
        std::cout << "Example: ./labosm osmtogeojson input.osm.pbf output.geojson" << std::endl;
        std::cout << "Example: ./labosm osmtofmi input.osm.pbf output.fmi\n" << std::endl;
        std::cout << "Example: ./labosm simpleserver input.fmi" << std::endl;
        std::cout << "Simple Server only supports dijkstra for routing\n" << std::endl;
        std::cout << "Example: ./labosm advancedserver input.fmi" << std::endl;
        std::cerr << "Advanced Server supports CH and Hub Labeling for routing" << std::endl;
        return 1;
    }

    if (argv[1] == std::string("simpleserver")) {
        if (argc != 3) {
            std::cerr << "Usage: ./labosm simpleserver input.fmi" << std::endl;
            return 1;
        }
        simpleServer(argv[2]);
        return 0;
    } else if (argv[1] == std::string("advancedserver")) {
        if (argc != 3) {
            std::cerr << "Usage: ./labosm advancedserver input.fmi" << std::endl;
            return 1;
        }
        advancedServer(argv[2]);
        return 0;
    }

    labosm::Graph g("../stgtregbz.fmi", true, 8, labosm::Heuristic::MIXED);

    // int dist = labosm::Graph::dijkstraQuery(g.getGraph(), 377371, 754742);
    // std::cout << "Distance: " << dist << std::endl;

    {
        labosm::Graph test("../stgtregbz.fmi", false, 1, labosm::Heuristic::IN_OUT);
        labosm::DijkstraQueryData data(test.getNumNodes());
        data.m_start = 377371;
        data.m_end = 754742;
        auto start = std::chrono::steady_clock::now();
        test.dijkstraQuery(data);
        auto end = std::chrono::steady_clock::now();
        std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
                  << std::endl;
        std::cout << "Distance: " << data.m_distance << std::endl;
    }

    {
        labosm::QueryData bd_data(g.getNumNodes());
        bd_data.m_start = 377371;
        bd_data.m_end = 754742;
        auto start = std::chrono::steady_clock::now();
        g.contractionHierarchyQuery(bd_data);
        auto end = std::chrono::steady_clock::now();
        std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
                  << std::endl;
        std::cout << "Distance: " << bd_data.m_distance << std::endl;
        std::cout << bd_data.m_meeting_node << std::endl;
        // g.bidirectionalDijkstraGetPath(bd_data);
    }

    /*
    g.createHubLabelsWithoutIS();
    {
        labosm::QueryData bd_data(g.getNumNodes());
        bd_data.m_start = 377371;
        bd_data.m_end = 754742;
        g.hubLabelQuery(bd_data);
        std::cout << "Distance: " << bd_data.m_distance << std::endl;
        std::cout << bd_data.m_meeting_node << std::endl;
    }


    std::cout << "Average Label size: " << g.averageLabelSize() << std::endl;
    std::cout << "Max Label size: " << g.maxLabelSize() << std::endl;
    */

    httplib::Server svr;
    // first the static content: website etc.
    svr.set_mount_point("/", "../static");

    // TODO: api for the graph
    // TODO: should work with JSON
    svr.Get(R"(/api/graph)", [&](const httplib::Request& req, httplib::Response& res) {
        res.set_content(R"({"status": "ok"})", "application/json");
    });

    // TODO: More

    std::cout << "Server started on port 8080" << std::endl;
    svr.listen("0.0.0.0", 8080);

    return 0;
}