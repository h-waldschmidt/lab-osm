#include "graph_creator.h"

#include <omp.h>

#include <fstream>
#include <iostream>
#include <random>

#define STB_IMAGE_IMPLEMENTATION
#include "../third-party/stb_image.h"

namespace labosm {
void CoastlineHandler::way(const osmium::Way& way, NodeMap& coastline_nodes, WayList& coastline_ways) {
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

void GraphCreator::generatePoints(const std::string& coastlines, int num_points, const std::string& output_file_prefix,
                                  bool image_based_filtering, const std::string& image_path) {
    auto start_reading = std::chrono::steady_clock::now();

    osmium::io::Reader reader{coastlines, osmium::osm_entity_bits::node | osmium::osm_entity_bits::way};

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

    write_geojson(output_file_prefix + "_coastlines.geojson", m_coastline_nodes, m_coastline_ways);

    auto points = generatePointsOnSphere(num_points);

    auto start_filtering = std::chrono::steady_clock::now();

    std::vector<std::pair<double, double>> filtered_points;
    if (image_based_filtering) {
        filtered_points = filterOutsideWaterImageBased(points, image_path);
    } else {
        filtered_points = filterOutsideWater(points);
    }

    auto end_filtering = std::chrono::steady_clock::now();
    std::cout << "Filtering time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_filtering - start_filtering).count()
              << " ms\n";
    std::cout << "Filtered points: " << filtered_points.size() << '\n';
    std::cout << "Original points: " << points.size() << '\n';
    std::cout << "Filtered points ratio: " << (double)filtered_points.size() / points.size() << '\n';

    std::ofstream out(output_file_prefix + "_filtered_points.geojson");
    out << R"({"type": "FeatureCollection", "features": [)";
    for (size_t i = 0; i < filtered_points.size(); ++i) {
        if (i > 0) out << ",";
        // point
        out << R"({"type": "Feature","geometry":{"type": "Point","coordinates":[)";
        out << filtered_points[i].first << "," << filtered_points[i].second << "]},";
        out << R"("properties":{}})";
    }
    out << "]}" << '\n';
    out.close();
    std::cout << "Filtered points written to " << output_file_prefix + "_filtered_points.geojson" << '\n';
}

void GraphCreator::generateGraph(const std::string& points_file, const std::string& output_file_path) {
    auto points = readPointsGeoJSON(points_file);
    std::cout << "Read " << points.size() << " points from geojson file\n";
    auto start = std::chrono::steady_clock::now();
    auto graph = createGraph(points);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Graph creation took: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms\n";
    printGraphToFMI(points, graph, output_file_path);
}

// https://mathworld.wolfram.com/SpherePointPicking.html
std::vector<std::pair<double, double>> GraphCreator::generatePointsOnSphere(int num_points) {
    std::vector<std::pair<double, double>> points(num_points);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    for (int i = 0; i < num_points; ++i) {
        double u = dis(gen);
        double v = dis(gen);

        double theta = 2 * M_PI * u;
        double phi = acos((2 * v) - 1);

        double lat = M_PI / 2 - phi;
        lat = lat * 180 / M_PI;
        double lon = theta * 180 / M_PI;
        lon = lon - 180;

        points[i] = {lon, lat};
    }

    return points;
}

void GraphCreator::write_geojson(const std::string& output_file, const NodeMap& nodes, const WayList& ways) {
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
            uint64_t last_id = it1->second.back();
            if (last_id == it1->first) continue;

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

// https://github.com/lcx366/SphericalPolygon/blob/master/sphericalpolygon/inside_polygon.py
std::vector<std::pair<double, double>> GraphCreator::filterOutsideWater(
    const std::vector<std::pair<double, double>>& points) {
    // generate bouding boxes for the ways
    std::vector<labosm::BoudingBox> boxes;
    boxes.reserve(m_coastline_ways.size());
    for (const auto& [first, way] : m_coastline_ways) {
        double min_lat = std::numeric_limits<double>::max();
        double max_lat = std::numeric_limits<double>::lowest();
        double min_lon = std::numeric_limits<double>::max();
        double max_lon = std::numeric_limits<double>::lowest();

        for (const auto& id : way) {
            const auto& [lon, lat] = m_coastline_nodes[id];
            min_lat = std::min(min_lat, lat);
            max_lat = std::max(max_lat, lat);
            min_lon = std::min(min_lon, lon);
            max_lon = std::max(max_lon, lon);
        }
        boxes.push_back(labosm::BoudingBox(min_lat, max_lat, min_lon, max_lon));
    }

    // convert points to 3D vectors
    std::vector<labosm::Vec3> vec_points(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& [lon, lat] = points[i];
        vec_points[i] = labosm::Vec3(lat, lon).normalize();
    }

    int max_way_size = 0;

    // convert ways to 3D vectors
    std::vector<std::vector<labosm::Vec3>> vec_ways;
    vec_ways.reserve(m_coastline_ways.size());
    for (const auto& [first, way] : m_coastline_ways) {
        vec_ways.push_back(std::vector<labosm::Vec3>(way.size()));
        auto& cache = vec_ways.back();

        max_way_size = std::max(max_way_size, (int)way.size());

        for (int j = 0; j < way.size(); ++j) {
            const auto& [lon, lat] = m_coastline_nodes[way[j]];
            labosm::Vec3 vec = labosm::Vec3(lat, lon).normalize();
            cache[j] = vec;
        }
    }

    int num_threads = 16;
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

// Projection: Equirectangular (aka Plate Carrée, EPSG:4326)
std::pair<int, int> GraphCreator::latlon_to_pixel(double lat, double lon, int img_width, int img_height) {
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

std::vector<std::pair<double, double>> GraphCreator::filterOutsideWaterImageBased(
    const std::vector<std::pair<double, double>>& points, const std::string& image_path) {
    int width, height, channels;
    unsigned char* data = stbi_load(image_path.c_str(), &width, &height, &channels, 0);
    if (!data) {
        std::cerr << "Failed to load image: " << image_path << "\n";
        return points;
    }
    std::vector<std::pair<double, double>> filtered_points;

    for (int i = 0; i < points.size(); ++i) {
        const auto& [lon, lat] = points[i];

        auto [px, py] = latlon_to_pixel(lat, lon, width, height);

        // skip if white
        if (data[py * width * channels + px * channels] == 255 &&
            data[py * width * channels + px * channels + 1] == 255 &&
            data[py * width * channels + px * channels + 2] == 255) {
            continue;
        }

        filtered_points.push_back(points[i]);
    }

    stbi_image_free(data);
    return filtered_points;
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

    /* ---------- 2. build point cloud in Cartesian coords ------------ */
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

    /* ---------- 3. build KD-tree (single threaded, few ms) ---------- */
    KDTree index(3, cloud, {10 /*max leaf*/});
    index.buildIndex();

    /* ---------- 4. radius search in parallel ----------------------- */
    std::vector<std::vector<labosm::Edge>> graph(N);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < N; ++i) {
        // 4.1  query neighbours inside 30 km chord
        const double q[3] = {cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z};
        std::vector<std::pair<size_t, double>> hits;
        nanoflann::SearchParams p;
        p.sorted = false;  // unsorted = faster
        index.radiusSearch(q, chord2, hits, p);

        // 4.2  keep closest in each quadrant
        int bestIdx[4] = {-1, -1, -1, -1};  // SW,SE,NW,NE
        int bestDist[4] = {INT_MAX, INT_MAX, INT_MAX, INT_MAX};

        const auto [lon1, lat1] = pointsDeg[i];
        const auto [lonRad, latRad] = pointsRad[i];
        for (auto& h : hits) {
            size_t j = h.first;
            if (j == i) continue;  // skip self

            int dist = greatCircleDistance(latRad, lonRad, pointsRad[j].second, pointsRad[j].first);
            if (dist > MAX_EDGE) continue;  // false-positive filter

            // quadrant test  (add constants to avoid  –180 / +180 wrap)
            double lat2a = pointsDeg[j].second + 90.0;
            double lon2a = pointsDeg[j].first + 180.0;
            double lat1a = lat1 + 90.0;
            double lon1a = lon1 + 180.0;

            bool south = lat2a < lat1a;
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

        if ((i & 0x3FFF) == 0)  // every ~16 k pts
            std::cout << "thread " << omp_get_thread_num() << " processed " << i << " vertices\n";
    }

    // backward edges
    for (size_t i = 0; i < N; ++i) {
        for (const auto& e : graph[i])  // forward edge (i → j)
        {
            const int j = e.m_target;

            bool alreadyThere = false;
            for (const auto& rev : graph[j])
                if (rev.m_target == static_cast<int>(i)) {  // reverse exists?
                    alreadyThere = true;
                    break;
                }
            if (alreadyThere) continue;  // skip duplicate

            graph[j].emplace_back(static_cast<int>(i), e.m_cost);
        }
    }

    return graph;
}

void GraphCreator::printGraphToFMI(const std::vector<std::pair<double, double>>& points,
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