#include "graph_creator.h"

#include <omp.h>

#include <algorithm>  // For std::min, std::max
#include <atomic>
#include <chrono>
#include <cmath>  // For M_PI, acos, sin, cos, atan2, fabs
#include <fstream>
#include <iostream>
#include <limits>  // For std::numeric_limits
#include <random>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include "../third-party/stb_image.h"

namespace labosm {

void GraphCreator::generatePointsAndFilter(const std::string& coastlines, int num_points,
                                           const std::string& output_file_prefix, bool image_based_filtering,
                                           const std::string& image_path) {
    if (!image_based_filtering) {
        auto start_reading = std::chrono::steady_clock::now();

        // Extract coastlines from OSM data using osmium
        osmium::io::Reader reader{coastlines, osmium::osm_entity_bits::node | osmium::osm_entity_bits::way};

        using index_type = osmium::index::map::SparseMemArray<osmium::unsigned_object_id_type, osmium::Location>;
        index_type index;
        osmium::handler::NodeLocationsForWays<index_type> location_handler{index};
        location_handler.ignore_errors();

        CoastlineHandler handler(&m_coastline_nodes, &m_coastline_ways);

        osmium::apply(reader, location_handler, handler);
        reader.close();
        auto end_reading = std::chrono::steady_clock::now();

        std::cout << "Reading time: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end_reading - start_reading).count()
                  << " ms\n";

        auto start_merging = std::chrono::steady_clock::now();
        merge_ways();
        auto end_merging = std::chrono::steady_clock::now();
        std::cout << "Merging time: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end_merging - start_merging).count()
                  << " ms\n";

        writeCoastlinesToGeojson(output_file_prefix + "_coastlines.geojson", m_coastline_nodes, m_coastline_ways);
    }

    std::vector<std::pair<double, double>> filtered_points;
    if (num_points <= 0) {
        std::cout << "Requested 0 or negative points. No points will be generated.\n";
    } else {
        filtered_points.reserve(num_points);
    }

    std::atomic<int> collected_points_count = 0;

    auto start_gen_and_filter = std::chrono::steady_clock::now();

    // Precomputation for polygon-based filtering
    std::vector<labosm::BoudingBox> poly_boxes;
    std::vector<std::vector<labosm::Vec3>> poly_vec_ways;

    if (!image_based_filtering) {
        poly_boxes.reserve(m_coastline_ways.size());
        for (const auto& way_entry : m_coastline_ways) {
            const auto& way_node_ids = way_entry.second;
            double min_lat = std::numeric_limits<double>::max();
            double max_lat = std::numeric_limits<double>::lowest();
            double min_lon = std::numeric_limits<double>::max();
            double max_lon = std::numeric_limits<double>::lowest();

            for (const auto& node_id : way_node_ids) {
                auto it = m_coastline_nodes.find(node_id);
                if (it != m_coastline_nodes.end()) {
                    const auto& [lon, lat] = it->second;
                    min_lat = std::min(min_lat, lat);
                    max_lat = std::max(max_lat, lat);
                    min_lon = std::min(min_lon, lon);
                    max_lon = std::max(max_lon, lon);
                }
            }
            poly_boxes.emplace_back(min_lat, max_lat, min_lon, max_lon);
        }

        poly_vec_ways.reserve(m_coastline_ways.size());
        for (const auto& way_entry : m_coastline_ways) {
            const auto& way_node_ids = way_entry.second;
            std::vector<labosm::Vec3> current_way_vecs;
            current_way_vecs.reserve(way_node_ids.size());
            for (const auto& node_id : way_node_ids) {
                auto it = m_coastline_nodes.find(node_id);
                if (it != m_coastline_nodes.end()) {
                    const auto& [lon, lat] = it->second;
                    current_way_vecs.push_back(labosm::Vec3(lat, lon).normalize());
                }
            }
            poly_vec_ways.push_back(current_way_vecs);
        }
    }

    // Image data loading for image-based filtering
    int img_width = 0, img_height = 0, img_channels = 0;
    unsigned char* img_data = nullptr;
    if (image_based_filtering) {
        img_data = stbi_load(image_path.c_str(), &img_width, &img_height, &img_channels, 0);
        if (!img_data) {
            std::cerr << "Failed to load image: " << image_path << "\n";
            // Exit or throw, as image-based filtering cannot proceed
            return;
        }
    }

#pragma omp parallel
    {
        std::mt19937 gen(std::random_device{}() + omp_get_thread_num());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        const double deg2rad_const = M_PI / 180.0;
        const double two_pi_const = 2.0 * M_PI;

        while (collected_points_count.load(std::memory_order_relaxed) < num_points) {
            double u = dis(gen);
            double v = dis(gen);
            double theta = two_pi_const * u;
            double phi = acos((2.0 * v) - 1.0);
            double lat_deg = (M_PI / 2.0 - phi) * 180.0 / M_PI;
            double lon_deg = (theta * 180.0 / M_PI) - 180.0;
            std::pair<double, double> current_point = {lon_deg, lat_deg};
            bool is_in_water = false;

            if (image_based_filtering) {
                is_in_water = isPointInWaterImageBased(lat_deg, lon_deg, img_data, img_width, img_height, img_channels);
            } else {  // Polygon-based filtering
                is_in_water = isPointInWaterPolygonBased(lat_deg, lon_deg, poly_vec_ways, poly_boxes);
            }

            if (is_in_water) {
#pragma omp critical
                {
                    if (filtered_points.size() < static_cast<size_t>(num_points)) {
                        filtered_points.push_back(current_point);
                        int current_total = collected_points_count.fetch_add(1, std::memory_order_relaxed) + 1;
                        if (current_total % 100000 == 0 || current_total == num_points) {
                            std::cout << "Collected " << current_total << " / " << num_points << " points...\n";
                        }
                    }
                }
            }
        }
    }

    if (image_based_filtering && img_data) {
        stbi_image_free(img_data);
    }

    auto end_gen_and_filter = std::chrono::steady_clock::now();
    std::cout
        << "Generation and Filtering time: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_gen_and_filter - start_gen_and_filter).count()
        << " ms\n";
    std::cout << "Target points: " << num_points << "\n";
    std::cout << "Collected points: " << filtered_points.size() << "\n";

    std::ofstream out(output_file_prefix + "_filtered_points.geojson");
    out << R"({"type": "FeatureCollection", "features": [)";
    for (size_t i = 0; i < filtered_points.size(); ++i) {
        if (i > 0) out << ",";
        out << R"({"type": "Feature","geometry":{"type": "Point","coordinates":[)";
        out << filtered_points[i].first << "," << filtered_points[i].second << "]},";
        out << R"("properties":{}})";
    }
    out << "]}" << '\n';
    out.close();
    std::cout << "Filtered points written to " << output_file_prefix + "_filtered_points.geojson" << '\n';
}

bool GraphCreator::isPointInWaterImageBased(double lat_deg, double lon_deg, const unsigned char* img_data,
                                            int img_width, int img_height, int img_channels) const {
    auto [px, py] = latlon_to_pixel(lat_deg, lon_deg, img_width, img_height);
    // Assuming white (255,255,255) is land
    if (img_channels >= 3 &&
        !(img_data[static_cast<size_t>(py) * img_width * img_channels + static_cast<size_t>(px) * img_channels] ==
              255 &&
          img_data[static_cast<size_t>(py) * img_width * img_channels + static_cast<size_t>(px) * img_channels + 1] ==
              255 &&
          img_data[static_cast<size_t>(py) * img_width * img_channels + static_cast<size_t>(px) * img_channels + 2] ==
              255)) {
        return true;
    } else if (img_channels == 1 && /* simple grayscale check, adapt if needed */
               !(img_data[static_cast<size_t>(py) * img_width + static_cast<size_t>(px)] == 255)) {
        return true;
    }
    // Add more conditions if other channel counts or color interpretations are needed
    return false;
}

bool GraphCreator::isPointInWaterPolygonBased(double lat_deg, double lon_deg,
                                              const std::vector<std::vector<labosm::Vec3>>& poly_vec_ways,
                                              const std::vector<labosm::BoudingBox>& poly_boxes) const {
    if (lat_deg <= -85.0) {  // South pole land (Antarctica)
        return false;        // Not in water
    }

    const double deg2rad_const = M_PI / 180.0;
    const double two_pi_const = 2.0 * M_PI;
    // labosm::Vec3 point_vec3 = labosm::Vec3(lat_deg, lon_deg).normalize(); // Not directly needed here, rotation
    // handles it
    bool inside_polygon = false;

    double r_theta_z = -lon_deg * deg2rad_const;
    double r_theta_y = (lat_deg - 90.0) * deg2rad_const;
    double cz = std::cos(r_theta_z);
    double sz = std::sin(r_theta_z);
    double cy = std::cos(r_theta_y);
    double sy = std::sin(r_theta_y);
    double R[3][3] = {{cy * cz, -cy * sz, sy}, {sz, cz, 0.0}, {-sy * cz, sy * sz, cy}};

    for (size_t j = 0; j < poly_vec_ways.size(); ++j) {
        // if (collected_points_count.load(std::memory_order_relaxed) >= num_points) break; // This check belongs in the
        // calling loop

        const auto& current_poly_vec_way = poly_vec_ways[j];
        if (current_poly_vec_way.empty()) continue;

        const auto& box = poly_boxes[j];
        if (!box.contains(lat_deg, lon_deg)) {
            continue;
        }

        double sum_angle = 0.0;
        double previous_lon_transformed = 0.0;

        const auto& first_node_vec = current_poly_vec_way[0];
        double x_rotated = R[0][0] * first_node_vec.x + R[0][1] * first_node_vec.y + R[0][2] * first_node_vec.z;
        double y_rotated = R[1][0] * first_node_vec.x + R[1][1] * first_node_vec.y + R[1][2] * first_node_vec.z;
        previous_lon_transformed = atan2(y_rotated, x_rotated);

        for (size_t k = 1; k < current_poly_vec_way.size(); ++k) {
            const auto& current_node_vec = current_poly_vec_way[k];
            x_rotated = R[0][0] * current_node_vec.x + R[0][1] * current_node_vec.y + R[0][2] * current_node_vec.z;
            y_rotated = R[1][0] * current_node_vec.x + R[1][1] * current_node_vec.y + R[1][2] * current_node_vec.z;

            double cur_lon_transformed = atan2(y_rotated, x_rotated);
            double angle = cur_lon_transformed - previous_lon_transformed;
            angle += (angle < -M_PI) * two_pi_const;
            angle -= (angle > M_PI) * two_pi_const;
            sum_angle += angle;
            previous_lon_transformed = cur_lon_transformed;
        }
        if (fabs(sum_angle - two_pi_const) < 1e-4 || fabs(sum_angle + two_pi_const) < 1e-4) {
            inside_polygon = true;
            break;
        }
    }
    return !inside_polygon;  // If not inside any polygon, it's in water
}

void GraphCreator::generateGraph(const std::string& points_file, const std::string& output_file_path) {
    auto points = readPointsGeoJSON(points_file);
    std::cout << "Read " << points.size() << " points from geojson file\n";
    auto start = std::chrono::steady_clock::now();
    auto graph = createGraph(points);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Graph creation took: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms\n";
    // At the end, randomize node order using DFS
    randomizeNodeOrderDFS(points, graph);
    writeGraphToFMI(points, graph, output_file_path);
}

void GraphCreator::randomizeNodeOrderDFS(std::vector<std::pair<double, double>>& points,
                                         std::vector<std::vector<labosm::Edge>>& graph) {
    size_t num_nodes = points.size();
    if (num_nodes == 0) return;
    std::vector<bool> visited(num_nodes, false);
    std::vector<int> order;
    order.reserve(num_nodes);
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, num_nodes - 1);
    size_t start_node = dist(rng);

    // Iterative DFS using explicit stack
    std::vector<size_t> stack;
    stack.push_back(start_node);
    while (!stack.empty()) {
        size_t node = stack.back();
        stack.pop_back();
        if (visited[node]) continue;
        visited[node] = true;
        order.push_back(static_cast<int>(node));
        // Push neighbors in reverse order for similar order as recursive
        const auto& edges = graph[node];
        for (auto it = edges.rbegin(); it != edges.rend(); ++it) {
            if (!visited[it->m_target]) {
                stack.push_back(it->m_target);
            }
        }
    }
    // If graph is disconnected, visit all components
    for (size_t i = 0; i < num_nodes; ++i) {
        if (!visited[i]) {
            stack.push_back(i);
            while (!stack.empty()) {
                size_t node = stack.back();
                stack.pop_back();
                if (visited[node]) continue;
                visited[node] = true;
                order.push_back(static_cast<int>(node));
                const auto& edges = graph[node];
                for (auto it = edges.rbegin(); it != edges.rend(); ++it) {
                    if (!visited[it->m_target]) {
                        stack.push_back(it->m_target);
                    }
                }
            }
        }
    }

    // Build old->new index map
    std::vector<int> old_to_new(num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) old_to_new[order[i]] = static_cast<int>(i);

    // Reorder points
    std::vector<std::pair<double, double>> new_points(num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) new_points[i] = points[order[i]];
    points = std::move(new_points);

    // Reorder graph and remap edges
    std::vector<std::vector<labosm::Edge>> new_graph(num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        for (const auto& edge : graph[order[i]]) {
            labosm::Edge new_edge = edge;
            new_edge.m_target = old_to_new[edge.m_target];
            new_graph[i].push_back(new_edge);
        }
    }
    graph = std::move(new_graph);
}

void GraphCreator::writeCoastlinesToGeojson(const std::string& output_file, const NodeMap& nodes, const WayList& ways) {
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

void GraphCreator::merge_ways() {
    size_t old_size = m_coastline_ways.size();
    while (true) {
        for (auto it1 = m_coastline_ways.begin(); it1 != m_coastline_ways.end(); ++it1) {
            // take the last node of the way
            uint64_t last_id = it1->second.back();
            if (last_id == it1->first) continue;

            // find the way that has the last node as first node and merge them
            auto it2 = m_coastline_ways.find(last_id);
            if (it2 != m_coastline_ways.end()) {
                it1->second.insert(it1->second.end(), it2->second.begin(), it2->second.end());
                m_coastline_ways.erase(it2);
            }
        }

        if (m_coastline_ways.size() == old_size || m_coastline_ways.size() == 1) {
            break;
        }

        old_size = m_coastline_ways.size();
    }
}

// Projection: Equirectangular (aka Plate Carrée)
std::pair<int, int> GraphCreator::latlon_to_pixel(double lat, double lon, int img_width, int img_height) const {
    // Normalize longitude from [−180, 180] to [0, 1]
    double x = (lon + 180.0) / 360.0;

    // Normalize latitude from [90, -90] to [0, 1] (top to bottom)
    double y = (90.0 - lat) / 180.0;

    // Convert to pixel coordinates
    int px = static_cast<int>(x * img_width);
    int py = static_cast<int>(y * img_height);

    // Clamp to image bounds just in case
    px = std::max(0, std::min(px, img_width - 1));
    py = std::max(0, std::min(py, img_height - 1));

    return {px, py};
}

std::vector<std::pair<double, double>> GraphCreator::readPointsGeoJSON(const std::string& filename) {
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

std::vector<std::vector<labosm::Edge>> GraphCreator::createGraph(std::vector<std::pair<double, double>>& pointsDeg) {
    constexpr int MAX_EDGE = 30'000;                                  // 30 km
    const double thetaMax = MAX_EDGE / R_EARTH;                       // rad
    const double chord2 = 4 * std::pow(std::sin(thetaMax * 0.5), 2);  // unit-sphere²

    // build point cloud with 3D coordinates
    const size_t N = pointsDeg.size();
    std::vector<std::pair<double, double>> pointsRad(N);
    PointCloud cloud;
    cloud.pts.resize(N);

#pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        auto [lonDeg, latDeg] = pointsDeg[i];
        double lat = latDeg * M_PI / 180.0, lon = lonDeg * M_PI / 180.0;
        pointsRad[i] = {lon, lat};

        double clat = std::cos(lat);
        cloud.pts[i] = {std::cos(lon) * clat, std::sin(lon) * clat, std::sin(lat)};
    }

    // build KDTree index
    KDTree index(3, cloud, {10 /*max leaf*/});
    index.buildIndex();

    // radius search
    std::vector<std::vector<labosm::Edge>> graph(N);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < N; ++i) {
        // query neighbours inside 30 km chord
        const double q[3] = {cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z};
        std::vector<std::pair<size_t, double>> hits;
        nanoflann::SearchParams p;
        p.sorted = false;  // unsorted = faster
        index.radiusSearch(q, chord2, hits, p);

        //  keep closest in each quadrant
        int bestIdx[4] = {-1, -1, -1, -1};  // SW,SE,NW,NE
        int bestDist[4] = {INT_MAX, INT_MAX, INT_MAX, INT_MAX};

        const auto [lon1, lat1] = pointsDeg[i];
        const auto [lonRad, latRad] = pointsRad[i];
        for (auto& h : hits) {
            size_t j = h.first;
            if (j == i) continue;  // skip self

            // false-positive filter
            int dist = greatCircleDistance(latRad, lonRad, pointsRad[j].second, pointsRad[j].first);
            if (dist > MAX_EDGE) continue;

            // quadrant test  (add constants to avoid  –180 / +180 wrap)
            double lat2a = pointsDeg[j].second + 90.0;
            double lon2a = pointsDeg[j].first + 180.0;
            double lat1a = lat1 + 90.0;
            double lon1a = lon1 + 180.0;

            // south test is easy
            bool south = lat2a < lat1a;
            // west test either by longitude or by wrapping
            bool west = (lon2a < lon1a) || (lon1a < 60.0 && lon2a > 300.0);
            int quadrant = south ? (west ? 0 : 1) : (west ? 2 : 3);

            if (dist < bestDist[quadrant]) {
                bestDist[quadrant] = dist;
                bestIdx[quadrant] = static_cast<int>(j);
            }
        }

        for (int q = 0; q < 4; ++q)
            if (bestIdx[q] != -1) {
                int j = bestIdx[q], d = bestDist[q];
                graph[i].emplace_back(j, d);
            }

        if (i % 10000 == 0) std::cout << "thread " << omp_get_thread_num() << " processed " << i << " vertices\n";
    }

    // backward edges
    for (size_t i = 0; i < N; ++i) {
        for (const auto& e : graph[i])  // forward edge (i → j)
        {
            const int j = e.m_target;

            bool alreadyThere = false;
            for (const auto& rev : graph[j])
                if (rev.m_target == static_cast<int>(i)) {
                    alreadyThere = true;
                    break;
                }
            if (alreadyThere) continue;  // skip duplicate

            graph[j].emplace_back(static_cast<int>(i), e.m_cost);
        }
    }

    return graph;
}

void GraphCreator::writeGraphToFMI(const std::vector<std::pair<double, double>>& points,
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
}  // namespace labosm