#include <omp.h>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
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
// #include "helper.h"

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

std::vector<std::pair<double, double>> generatePointsOnSphere(int num_points) {
    std::vector<std::pair<double, double>> points(num_points);

    // fibonacci sphere algorithm
    // https://openprocessing.org/sketch/41142

    double phi = (sqrt(5) + 1) / 2 - 1;  // golden ratio
    double theta = 2 * M_PI * phi;       // golden angle

    for (int i = 0; i < num_points; ++i) {
        // long must be between -pi and pi
        // lat must be between -pi/2 and pi/2
        double lat = asin(-1 + (2 * i / (double)num_points));
        double lon = theta * i;
        lon /= 2 * M_PI;
        lon -= floor(lon);
        lon *= 2 * M_PI;
        if (lon > M_PI) lon -= 2 * M_PI;

        // transform to degrees
        lat = lat * 180 / M_PI;
        lon = lon * 180 / M_PI;

        std::cout << lat << " " << lon << std::endl;

        points[i] = {lon, lat};
    }

    return points;
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

// https://github.com/lcx366/SphericalPolygon/blob/master/sphericalpolygon/inside_polygon.py
std::vector<std::pair<double, double>> filterOutsideWater(const std::vector<std::pair<double, double>>& points) {
    // generate bouding boxes for the ways
    std::vector<BoudingBox> boxes;
    boxes.reserve(coastline_ways.size());
    for (const auto& [first, way] : coastline_ways) {
        double min_lat = std::numeric_limits<double>::max();
        double max_lat = std::numeric_limits<double>::lowest();
        double min_lon = std::numeric_limits<double>::max();
        double max_lon = std::numeric_limits<double>::lowest();

        for (const auto& id : way) {
            const auto& [lon, lat] = coastline_nodes[id];
            min_lat = std::min(min_lat, lat);
            max_lat = std::max(max_lat, lat);
            min_lon = std::min(min_lon, lon);
            max_lon = std::max(max_lon, lon);
        }
        boxes.push_back(BoudingBox(min_lat, max_lat, min_lon, max_lon));
    }

    // convert points to 3D vectors
    std::vector<Vec3> vec_points(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& [lon, lat] = points[i];
        vec_points[i] = Vec3(lat, lon).normalize();
    }

    int max_way_size = 0;

    // convert ways to 3D vectors
    std::vector<std::vector<Vec3>> vec_ways;
    vec_ways.reserve(coastline_ways.size());
    for (const auto& [first, way] : coastline_ways) {
        vec_ways.push_back(std::vector<Vec3>(way.size()));
        auto& cache = vec_ways.back();

        max_way_size = std::max(max_way_size, (int)way.size());

        for (int j = 0; j < way.size(); ++j) {
            const auto& [lon, lat] = coastline_nodes[way[j]];
            Vec3 vec = Vec3(lat, lon).normalize();
            cache[j] = vec;
        }
    }

    int num_threads = 4;
    omp_set_num_threads(num_threads);
    std::vector<std::vector<int>> inside_points_per_thread(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        inside_points_per_thread[i].reserve(points.size() / num_threads);
    }

    const double deg2rad = M_PI / 180.0;
    const double two_pi = 2 * M_PI;
    // check if points are inside one of the polygons (if yes they are not in the water and must be filtered out)
#pragma omp parallel for
    for (int i = 0; i < points.size(); ++i) {
        auto& inside_points = inside_points_per_thread[omp_get_thread_num()];

        const auto& [lon, lat] = points[i];

        // south pole edge cases that the bounding box test doesnt cover
        if (lat <= -85) {
            continue;
        }

        const auto& point = vec_points[i];
        bool inside = false;

        double theta_z = -lon * deg2rad;          // rotation around z-axis
        double theta_y = (lat - 90.0) * deg2rad;  // rotation around y-axis

        // Precompute sines and cosines.
        double cz = std::cos(theta_z);
        double sz = std::sin(theta_z);
        double cy = std::cos(theta_y);
        double sy = std::sin(theta_y);

        // Build the composite rotation matrix R = Ry * Rz.
        // The 3x3 matrix is computed as:
        // [ [cy*cz,    -cy*sz,    sy],
        //   [sz,       cz,        0],
        //   [-sy*cz,   sy*sz,     cy] ]
        // clang-format off
        double R[3][3] = {
            {cy * cz, -cy * sz, sy},
            {sz,      cz,       0.0},
            {-sy * cz, sy * sz, cy}
        };
        // clang-format on

        for (size_t j = 0; j < vec_ways.size(); ++j) {
            const auto& vec_way = vec_ways[j];
            const auto& box = boxes[j];
            // just early termination
            // can skip the way, because it is not in the bounding box

            if (!box.contains(lat, lon)) {
                continue;
            }

            double sum_angle = 0.0;
            double cur_lon_transformed = 0.0;
            double previous_lon_transformed = 0.0;

            // calculate the first iteration
            const auto& vec = vec_way[0];
            double x_rotated = R[0][0] * vec.x + R[0][1] * vec.y + R[0][2] * vec.z;
            double y_rotated = R[1][0] * vec.x + R[1][1] * vec.y + R[1][2] * vec.z;
            double z_rotated = R[2][0] * vec.x + R[2][1] * vec.y + R[2][2] * vec.z;
            previous_lon_transformed = atan2(y_rotated, x_rotated);

            // #pragma omp simd
            for (size_t k = 1; k < vec_way.size(); ++k) {
                const auto& vec = vec_way[k];
                // apply rotation
                x_rotated = R[0][0] * vec.x + R[0][1] * vec.y + R[0][2] * vec.z;
                y_rotated = R[1][0] * vec.x + R[1][1] * vec.y + R[1][2] * vec.z;
                z_rotated = R[2][0] * vec.x + R[2][1] * vec.y + R[2][2] * vec.z;

                cur_lon_transformed = atan2(y_rotated, x_rotated);

                double angle = cur_lon_transformed - previous_lon_transformed;
                angle += (angle < -M_PI) * two_pi;
                angle -= (angle > M_PI) * two_pi;
                sum_angle += angle;

                previous_lon_transformed = cur_lon_transformed;
            }

            if (fabs(sum_angle - two_pi) < 1e-4 || fabs(sum_angle + two_pi) < 1e-4) {
                inside = true;
                break;
            }
        }

        if (inside) {
            continue;
        }

        inside_points.push_back(i);

        if (i % 1000 == 0) {
            std::cout << "Thread " << omp_get_thread_num() << " processed " << inside_points.size() << " points\n";
        }
    }

    // collect all points from all threads
    std::vector<std::pair<double, double>> filtered_points;
    for (const auto& inside_points : inside_points_per_thread) {
        for (const auto& index : inside_points) {
            filtered_points.push_back(points[index]);
        }
    }

    return filtered_points;
}

std::vector<std::pair<double, double>> readPointsGeoJSON(const std::string& filename) {
    std::ifstream in(filename);
    std::string line;
    std::vector<std::pair<double, double>> points;

    // all the json is in one line
    std::getline(in, line);
    std::string::size_type pos = 0;
    while ((pos = line.find("\"coordinates\"", pos)) != std::string::npos) {
        pos = line.find("[", pos);
        if (pos == std::string::npos) break;

        std::string::size_type end_pos = line.find("]", pos);
        std::string coords = line.substr(pos + 1, end_pos - pos);

        std::stringstream ss(coords);
        double lon, lat;
        char comma;
        while (ss >> lon >> comma >> lat) {
            points.emplace_back(lon, lat);
            if (ss.peek() == ',') ss.ignore();
        }
        pos = end_pos + 1;
    }
    in.close();

    return points;
}

std::vector<std::vector<labosm::Edge>> createGraph(std::vector<std::pair<double, double>>& points) {
    const int max_edge_length = 300000000;  // 30 km
    std::vector<std::vector<labosm::Edge>> graph(points.size());

    for (int i = 0; i < points.size(); ++i) {
        const auto& [lon1, lat1] = points[i];
        // select closest points in these directions
        int south_west_index = -1;
        int south_west_distance = std::numeric_limits<int>::max();
        int south_east_index = -1;
        int south_east_distance = std::numeric_limits<int>::max();
        int north_west_index = -1;
        int north_west_distance = std::numeric_limits<int>::max();
        int north_east_index = -1;
        int north_east_distance = std::numeric_limits<int>::max();

        for (int j = 0; j < points.size(); ++j) {
            if (i == j) continue;

            const auto& [lon1, lat1] = points[i];
            const auto& [lon2, lat2] = points[j];

            double distance = labosm::Graph::greatCircleDistance(lat1, lon1, lat2, lon2);
            if (distance > max_edge_length) continue;

            // edge cases at -180, 180, -90, 90
            // add 180 and 90 to the lon and lat to reduce the edge cases
            double lat1_adjusted = lat1 + 90;
            double lon1_adjusted = lon1 + 180;
            double lat2_adjusted = lat2 + 90;
            double lon2_adjusted = lon2 + 180;

            if (lat2_adjusted < lat1_adjusted && lon2_adjusted > lon1_adjusted && lon1_adjusted < 180 &&
                lon2_adjusted > 180 && distance < south_west_distance) {
                south_west_distance = distance;
                south_west_index = j;
            } else if (lat2_adjusted < lat1_adjusted && lon2_adjusted < lon1_adjusted && lon1_adjusted > 180 &&
                       lon2_adjusted < 180 && distance < south_east_distance) {
                south_east_distance = distance;
                south_east_index = j;
            } else if (lat2_adjusted > lat1_adjusted && lon2_adjusted < lon1_adjusted && lon1_adjusted > 180 &&
                       lon2_adjusted < 180 && distance < north_west_distance) {
                north_west_distance = distance;
                north_west_index = j;
            } else if (lat2_adjusted > lat1_adjusted && lon2_adjusted > lon1_adjusted && lon1_adjusted < 180 &&
                       lon2_adjusted > 180 && distance < north_east_distance) {
                north_east_distance = distance;
                north_east_index = j;
            } else if (lat2_adjusted < lat1_adjusted && lon2_adjusted < lon1_adjusted &&
                       distance < south_west_distance) {
                south_west_distance = distance;
                south_west_index = j;
            } else if (lat2_adjusted < lat1_adjusted && lon2_adjusted > lon1_adjusted &&
                       distance < south_east_distance) {
                south_east_distance = distance;
                south_east_index = j;
            } else if (lat2_adjusted > lat1_adjusted && lon2_adjusted < lon1_adjusted &&
                       distance < north_west_distance) {
                north_west_distance = distance;
                north_west_index = j;
            } else if (lat2_adjusted > lat1_adjusted && lon2_adjusted > lon1_adjusted &&
                       distance < north_east_distance) {
                north_east_distance = distance;
                north_east_index = j;
            }
        }

        // check that none of the indices are the same
        if (south_west_index != -1) {
            graph[i].emplace_back(south_west_index, south_west_distance);
        }
        if (south_east_index != -1 && south_east_index != south_west_index) {
            graph[i].emplace_back(south_east_index, south_east_distance);
        }
        if (north_west_index != -1 && north_west_index != south_west_index && north_west_index != south_east_index) {
            graph[i].emplace_back(north_west_index, north_west_distance);
        }
        if (north_east_index != -1 && north_east_index != south_west_index && north_east_index != south_east_index &&
            north_east_index != north_west_index) {
            graph[i].emplace_back(north_east_index, north_east_distance);
        }
    }

    return graph;
}

void printGraphToFMI(const std::vector<std::pair<double, double>>& points,
                     const std::vector<std::vector<labosm::Edge>>& graph, const std::string& filename) {
    std::ofstream out(filename);
    out << "# Timestamp: " << std::time(nullptr) << '\n';
    out << "# Type: Coastlines \n\n";

    int num_edges = 0;
    for (const auto& edges : graph) {
        num_edges += edges.size();
    }

    out << points.size() << '\n';
    out << num_edges << '\n';

    for (int i = 0; i < points.size(); ++i) {
        const auto& [lon, lat] = points[i];
        out << i << " 0 " << lat << " " << lon << " 0" << '\n';
    }

    for (int i = 0; i < points.size(); ++i) {
        const auto& edges = graph[i];
        for (const auto& edge : edges) {
            out << i << " " << edge.m_target << " " << edge.m_cost << " 0 0\n";
        }
    }
}

/*
int main() {
    auto points = readPointsGeoJSON("filtered_points.geojson");
    std::cout << "Read " << points.size() << " points from geojson file\n";
    auto start = std::chrono::steady_clock::now();
    auto graph = createGraph(points);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Graph creation tooK: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms\n";
    printGraphToFMI(points, graph, "coastlines.fmi");
}
*/

int main(int argc, char* argv[]) {
    const char* input_filename = argv[1];

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

    // size_t max_size = 0;
    // std::vector<uint64_t> longest_way;
    // for (const auto& way : coastline_ways) {
    //     if (way.second.size() > max_size) {
    //         max_size = way.second.size();
    //         longest_way = way.second;
    //     }
    // }

    // coastline_ways.clear();
    // coastline_ways[longest_way.front()] = longest_way;

    int num_points = 1000;
    auto start_filtering = std::chrono::steady_clock::now();
    auto points = generatePointsOnSphere(num_points);

    std::ofstream out("points.geojson");
    out << R"({"type": "FeatureCollection", "features": [)";
    for (size_t i = 0; i < points.size(); ++i) {
        if (i > 0) out << ",";
        // point
        out << R"({"type": "Feature","geometry":{"type": "Point","coordinates":[)";
        out << points[i].first << "," << points[i].second << "]},";
        out << R"("properties":{}})";
    }
    out << "]}" << '\n';
    out.close();

    auto filtered_points = filterOutsideWater(points);
    auto end_filtering = std::chrono::steady_clock::now();
    std::cout << "Filtering time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_filtering - start_filtering).count()
              << " ms\n";

    std::ofstream out_filtered("filtered_points.geojson");
    out_filtered << R"({"type": "FeatureCollection", "features": [)";
    // write all into one LineString
    for (size_t i = 0; i < filtered_points.size(); ++i) {
        if (i > 0) out_filtered << ",";
        // point
        out_filtered << R"({"type": "Feature","geometry":{"type": "Point","coordinates":[)";
        out_filtered << filtered_points[i].first << "," << filtered_points[i].second << "]},";
        out_filtered << R"("properties":{}})";
    }
    out_filtered << "]}" << '\n';
    out_filtered.close();

    // write_geojson(output_filename, coastline_nodes, coastline_ways);

    return 0;
}

/*
int main() {
    // render as lines with point duplicated
    // make points much bigger
    std::ofstream out("points.geojson");
    out << R"({"type": "FeatureCollection", "features": [)";
    // write all into one LineString
    for (size_t i = 0; i < points.size(); ++i) {
        if (i > 0) out << ",";
        // point
        out << R"({"type": "Feature","geometry":{"type": "Point","coordinates":[)";
        out << points[i].first << "," << points[i].second << "]},";
        out << R"("properties":{}})";
    }
    out << "]}" << '\n';
    out.close();

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
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()
                      << std::endl;
            if (data.m_distance == std::numeric_limits<int>::max()) {
                res.status = 400;
                res.set_content(R"({"error": "No path found"})", "application/json");
                return;
            }

            start_time = std::chrono::steady_clock::now();
            g.dijkstraExtractPath(data);
            end_time = std::chrono::steady_clock::now();
            std::cout << "Dijkstra Path Extraction Time: "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()
                      << std::endl;

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
            res.set_content(R"({"error": "Missing start or end
parameter"})",
                            "application/json");
        }
    });

    std::cout << "Server started on port 8080" << std::endl;
    svr.listen("0.0.0.0", 8080);
}

void advancedServer(const std::string& fmi_file) {
    labosm::Graph g(fmi_file, true, 4, labosm::Heuristic::MIXED);
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
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()
                      << std::endl;
            std::cout << "Start: " << data.m_start << std::endl;
            std::cout << "End: " << data.m_end << std::endl;
            std::cout << "Distance: " << data.m_distance << std::endl;
            std::cout << "Meeting Node: " << data.m_meeting_node << std::endl;

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
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()
                      << std::endl;

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

            std::cout << "Distance: " << data.m_distance << std::endl;

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
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()
                      << std::endl;
            std::cout << "Start: " << data.m_start << std::endl;
            std::cout << "End: " << data.m_end << std::endl;
            std::cout << "Distance: " << data.m_distance << std::endl;
            std::cout << "Meeting Node: " << data.m_meeting_node << std::endl;

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
                      << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()
                      << std::endl;

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

    std::cout << "Server started on port 8080" << std::endl;
    svr.listen("0.0.0.0", 8080);
}
/*
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

    // labosm::Graph g("../stgtregbz.fmi", true, 8, labosm::Heuristic::MIXED);

    // // int dist = labosm::Graph::dijkstraQuery(g.getGraph(), 377371, 754742);
    // // std::cout << "Distance: " << dist << std::endl;

    // {
    //     labosm::Graph test("../stgtregbz.fmi", false, 1, labosm::Heuristic::IN_OUT);
    //     labosm::DijkstraQueryData data(test.getNumNodes());
    //     data.m_start = 377371;
    //     data.m_end = 754742;
    //     auto start = std::chrono::steady_clock::now();
    //     test.dijkstraQuery(data);
    //     auto end = std::chrono::steady_clock::now();
    //     std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
    //               << std::endl;
    //     std::cout << "Distance: " << data.m_distance << std::endl;
    // }

    // {
    //     labosm::QueryData bd_data(g.getNumNodes());
    //     bd_data.m_start = 377371;
    //     bd_data.m_end = 754742;
    //     auto start = std::chrono::steady_clock::now();
    //     g.contractionHierarchyQuery(bd_data);
    //     auto end = std::chrono::steady_clock::now();
    //     std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
    //               << std::endl;
    //     std::cout << "Distance: " << bd_data.m_distance << std::endl;
    //     std::cout << bd_data.m_meeting_node << std::endl;
    //     // g.bidirectionalDijkstraGetPath(bd_data);
    // }

    return 0;
}
*/