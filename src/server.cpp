#include "server.h"

#include "graph.h"
#include "httplib.h"

namespace labosm::server {

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
            auto query_start_time = std::chrono::steady_clock::now();
            g.dijkstraQuery(data);
            auto query_end_time = std::chrono::steady_clock::now();
            long long query_time_ms =
                std::chrono::duration_cast<std::chrono::microseconds>(query_end_time - query_start_time).count();
            std::cout << "Dijkstra Query Time: " << query_time_ms << "\n";
            std::cout << "Memory Usage: " << getMemoryUsage() << " KB\n";
            std::cout << "PQ Pops: " << data.num_pq_pops << "\n";
            if (data.m_distance == std::numeric_limits<int>::max()) {
                res.status = 400;
                res.set_content(R"({"error": "No path found"})", "application/json");
                return;
            }

            auto extraction_start_time = std::chrono::steady_clock::now();
            g.dijkstraExtractPath(data);
            auto extraction_end_time = std::chrono::steady_clock::now();
            long long extraction_time_ms =
                std::chrono::duration_cast<std::chrono::microseconds>(extraction_end_time - extraction_start_time)
                    .count();
            std::cout << "Dijkstra Path Extraction Time: " << extraction_time_ms << "\n";

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

            std::cout << "Distance: " << data.m_distance << "\n";

            res.set_content(R"({"distance": )" + std::to_string(data.m_distance) + R"(,"path": )" + path_str +
                                R"(,"pq_pops": )" + std::to_string(data.num_pq_pops) + R"(,"memory_usage": )" +
                                std::to_string(getMemoryUsage()) + R"(,"query_time": )" +
                                std::to_string(query_time_ms) + R"(,"extraction_time": )" +
                                std::to_string(extraction_time_ms) + R"(})",
                            "application/json");
        } else {
            res.status = 400;
            res.set_content(R"({"error": "Missing start or end parameter"})", "application/json");
        }
    });

    std::cout << "Server started on port 8080" << "\n";
    svr.listen("0.0.0.0", 8080);
}

void advancedServer(const std::string& fmi_file, bool enable_hub_labels) {
    labosm::Graph g(fmi_file, true, 16, labosm::Heuristic::MIXED);
    labosm::QueryData data(g.getNumNodes());

    if (enable_hub_labels) {
        g.createHubLabelsWithIS();
    }
    std::string enable_hub_labels_str = enable_hub_labels ? "true" : "false";

    httplib::Server svr;
    // first the static content: website etc.
    svr.set_mount_point("/", "static");

    svr.Get(R"(/api/)", [&](const httplib::Request& req, httplib::Response& res) {
        res.set_content(R"({"type": "complex", "hub_labels": )" + enable_hub_labels_str + R"(})", "application/json");
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
            auto query_start_time = std::chrono::steady_clock::now();
            g.contractionHierarchyQuery(data);
            auto query_end_time = std::chrono::steady_clock::now();
            long long query_time_ms =
                std::chrono::duration_cast<std::chrono::microseconds>(query_end_time - query_start_time).count();
            std::cout << "CH Query Time: " << query_time_ms << "\n";
            std::cout << "Memory Usage: " << getMemoryUsage() << " KB\n";
            std::cout << "PQ Pops: " << data.num_pq_pops << "\n";
            std::cout << "Start: " << data.m_start << "\n";
            std::cout << "End: " << data.m_end << "\n";
            std::cout << "Distance: " << data.m_distance << "\n";
            std::cout << "Meeting Node: " << data.m_meeting_node << "\n";

            if (data.m_meeting_node == -1) {
                res.status = 400;
                res.set_content(R"({"error": "No path found"})", "application/json");
                return;
            }

            auto extraction_start_time = std::chrono::steady_clock::now();
            // TODO: There still seems to be a bug in the path extraction
            // It gets stuck at one of the while loops
            // this occurs very rarely
            // I think this might be an issue due to multiple requests at the same time
            // So a mutex might be needed

            g.contractionHierarchyExtractPath(data);
            auto extraction_end_time = std::chrono::steady_clock::now();
            long long extraction_time_ms =
                std::chrono::duration_cast<std::chrono::microseconds>(extraction_end_time - extraction_start_time)
                    .count();
            std::cout << "CH Path Extraction Time: " << extraction_time_ms << "\n";

            res.status = 200;

            std::vector<std::pair<double, double>> coords;
            for (int node : data.m_shortest_path) {
                coords.push_back(g.getNodeCoords(node));
            }
            std::string path_str = "{\"type\": \"LineString\", \"coordinates\": [";
            for (size_t i = 0; i < coords.size(); ++i) {
                if (i > 0) path_str += ",";
                path_str += "[" + std::to_string(coords[i].first) + "," + std::to_string(coords[i].second) + "]";
            }
            path_str += "]}";

            std::cout << "Distance: " << data.m_distance << "\n";

            res.set_content(R"({"distance": )" + std::to_string(data.m_distance) + R"(,"path": )" + path_str +
                                R"(,"pq_pops": )" + std::to_string(data.num_pq_pops) + R"(,"memory_usage": )" +
                                std::to_string(getMemoryUsage()) + R"(,"query_time": )" +
                                std::to_string(query_time_ms) + R"(,"extraction_time": )" +
                                std::to_string(extraction_time_ms) + R"(})",
                            "application/json");
        } else {
            res.status = 400;
            res.set_content(R"({"error": "Missing start or end parameter"})", "application/json");
        }
    });

    // Hub Label Query
    svr.Get(R"(/api/hub-label)", [&](const httplib::Request& req, httplib::Response& res) {
        if (!enable_hub_labels) {
            res.status = 400;
            res.set_content(R"({"error": "Hub labeling is disabled on this server", "hub_labels": false})",
                            "application/json");
            return;
        }
        if (req.has_param("start") && req.has_param("end")) {
            int start = std::stoi(req.get_param_value("start"));
            int end = std::stoi(req.get_param_value("end"));
            data.m_start = start;
            data.m_end = end;
            auto query_start_time = std::chrono::steady_clock::now();
            auto hub_indices = g.hubLabelQuery(data);
            auto query_end_time = std::chrono::steady_clock::now();
            long long query_time_ms =
                std::chrono::duration_cast<std::chrono::microseconds>(query_end_time - query_start_time).count();
            std::cout << "Hub Label Query Time: " << query_time_ms << "\n";
            std::cout << "Memory Usage: " << getMemoryUsage() << " KB\n";
            std::cout << "PQ Pops: " << data.num_pq_pops << "\n";
            std::cout << "Start: " << data.m_start << "\n";
            std::cout << "End: " << data.m_end << "\n";
            std::cout << "Distance: " << data.m_distance << "\n";
            std::cout << "Meeting Node: " << data.m_meeting_node << "\n";

            if (data.m_meeting_node == -1) {
                res.status = 400;
                res.set_content(R"({"error": "No path found"})", "application/json");
                return;
            }

            auto extraction_start_time = std::chrono::steady_clock::now();
            // TODO: There still seems to be a bug in the path extraction
            // It gets stuck at one of the while loops
            // this occurs very rarely
            // I think this might be an issue due to multiple requests at the same time
            // So a mutex might be needed

            g.hubLabelExtractPath(data, hub_indices);
            auto extraction_end_time = std::chrono::steady_clock::now();
            long long extraction_time_ms =
                std::chrono::duration_cast<std::chrono::microseconds>(extraction_end_time - extraction_start_time)
                    .count();
            std::cout << "Hub Label Path Extraction Time: " << extraction_time_ms << "\n";

            res.status = 200;

            std::vector<std::pair<double, double>> coords;
            for (int node : data.m_shortest_path) {
                coords.push_back(g.getNodeCoords(node));
            }

            std::string path_str = "{\"type\": \"LineString\", \"coordinates\": [";
            for (size_t i = 0; i < coords.size(); ++i) {
                if (i > 0) path_str += ",";
                path_str += "[" + std::to_string(coords[i].first) + "," + std::to_string(coords[i].second) + "]";
            }
            path_str += "]}";

            res.set_content(R"({"distance": )" + std::to_string(data.m_distance) + R"(,"path": )" + path_str +
                                R"(,"pq_pops": )" + std::to_string(data.num_pq_pops) + R"(,"memory_usage": )" +
                                std::to_string(getMemoryUsage()) + R"(,"query_time": )" +
                                std::to_string(query_time_ms) + R"(,"extraction_time": )" +
                                std::to_string(extraction_time_ms) + R"(})",
                            "application/json");
        } else {
            res.status = 400;
            res.set_content(R"({"error": "Missing start or end parameter"})", "application/json");
        }
    });

    std::cout << "Server started on port 8080" << "\n";
    svr.listen("0.0.0.0", 8080);
}
}  // namespace labosm::server