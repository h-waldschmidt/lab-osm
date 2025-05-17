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
            res.set_content(R"({"error": "Missing lat or lon
parameter"})",
                            "application/json");
        }
    });

    // dijkstra query
    svr.Get(R"(/api/dijkstra)", [&](const httplib::Request& req, httplib::Response& res) {
        if (req.has_param("start") && req.has_param("end")) {
            int start = std::stoi(req.get_param_value("start"));
            int end = std::stoi(req.get_param_value("end"));
            data.m_start = start;
            data.m_end = end;
            auto start_time = std::chrono::steady_clock::now();
            g.dijkstraQuery(data);
            auto end_time = std::chrono::steady_clock::now();
            std::cout << "Dijkstra Query Time: "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() << "\n";
            if (data.m_distance == std::numeric_limits<int>::max()) {
                res.status = 400;
                res.set_content(R"({"error": "No path found"})", "application/json");
                return;
            }

            start_time = std::chrono::steady_clock::now();
            g.dijkstraExtractPath(data);
            end_time = std::chrono::steady_clock::now();
            std::cout << "Dijkstra Path Extraction Time: "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() << "\n";

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

            res.set_content(R"({"distance": )" + std::to_string(data.m_distance) + R"(,"path": )" + path_str + R"(})",
                            "application/json");
        } else {
            res.status = 400;
            res.set_content(R"({"error": "Missing start or end
parameter"})",
                            "application/json");
        }
    });

    std::cout << "Server started on port 8080" << "\n";
    svr.listen("0.0.0.0", 8080);
}

void advancedServer(const std::string& fmi_file) {
    labosm::Graph g(fmi_file, true, 12, labosm::Heuristic::MIXED);
    labosm::QueryData data(g.getNumNodes());

    g.createHubLabelsWithIS();

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
            res.set_content(R"({"error": "Missing lat or lon
parameter"})",
                            "application/json");
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
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() << "\n";
            std::cout << "Start: " << data.m_start << "\n";
            std::cout << "End: " << data.m_end << "\n";
            std::cout << "Distance: " << data.m_distance << "\n";
            std::cout << "Meeting Node: " << data.m_meeting_node << "\n";

            if (data.m_meeting_node == -1) {
                res.status = 400;
                res.set_content(R"({"error": "No path found"})", "application/json");
                return;
            }

            start_time = std::chrono::steady_clock::now();
            // TODO: There still seems to be a bug in the path extraction
            // It gets stuck at one of the while loops
            g.contractionHierarchyExtractPath(data);
            end_time = std::chrono::steady_clock::now();
            std::cout << "CH Path Extraction Time: "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() << "\n";

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

            res.set_content(R"({"distance": )" + std::to_string(data.m_distance) + R"(,"path": )" + path_str + R"(})",
                            "application/json");
        } else {
            res.status = 400;
            res.set_content(R"({"error": "Missing start or end parameter"})", "application/json");
        }
    });

    // Hub Label Query
    svr.Get(R"(/api/hub-label)", [&](const httplib::Request& req, httplib::Response& res) {
        if (req.has_param("start") && req.has_param("end")) {
            int start = std::stoi(req.get_param_value("start"));
            int end = std::stoi(req.get_param_value("end"));
            data.m_start = start;
            data.m_end = end;
            auto start_time = std::chrono::steady_clock::now();
            auto hub_indices = g.hubLabelQuery(data);
            auto end_time = std::chrono::steady_clock::now();
            std::cout << "Hub Label Query Time: "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() << "\n";
            std::cout << "Start: " << data.m_start << "\n";
            std::cout << "End: " << data.m_end << "\n";
            std::cout << "Distance: " << data.m_distance << "\n";
            std::cout << "Meeting Node: " << data.m_meeting_node << "\n";

            if (data.m_meeting_node == -1) {
                res.status = 400;
                res.set_content(R"({"error": "No path found"})", "application/json");
                return;
            }

            start_time = std::chrono::steady_clock::now();
            // TODO: There still seems to be a bug in the path extraction
            // It gets stuck at one of the while loops
            g.hubLabelExtractPath(data, hub_indices);
            end_time = std::chrono::steady_clock::now();
            std::cout << "Hub Label Path Extraction Time: "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() << "\n";

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

            res.set_content(R"({"distance": )" + std::to_string(data.m_distance) + R"(,"path": )" + path_str + R"(})",
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