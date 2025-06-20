#pragma once

#include <nanoflann.hpp>
#include <osmium/handler.hpp>
#include <osmium/handler/node_locations_for_ways.hpp>
#include <osmium/index/map/sparse_mem_array.hpp>
#include <osmium/io/any_input.hpp>
#include <osmium/osm/location.hpp>
#include <osmium/visitor.hpp>

#include "helper.h"

namespace labosm {
// Used for the extraction of the coastlines from the osm.pbf files
using NodeMap = std::unordered_map<uint64_t, std::pair<double, double>>;
using WayList = std::unordered_map<uint64_t, std::vector<uint64_t>>;

// KDTree for the graph creation
using KDTree =
    nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3, size_t>;

/**
 * @brief Handler for coastlines in OSM data.
 */
class CoastlineHandler : public osmium::handler::Handler {
   public:
    CoastlineHandler(NodeMap* nodes, WayList* ways) : coastline_nodes(nodes), coastline_ways(ways) {}
    ~CoastlineHandler() = default;

    /**
     * @brief Handle a way in the OSM data.
     * Save the coordinates of a node in a hashmap and the way in a vector.
     * The way is indexed by the first node of the way, to make merging easier.
     *
     * @param way The way to handle.
     */
    void way(const osmium::Way& way) {
        if (way.tags().has_tag("natural", "coastline")) {
            std::vector<uint64_t> way_nodes;
            for (const auto& node_ref : way.nodes()) {
                osmium::Location loc = node_ref.location();
                if (loc.valid()) {
                    uint64_t id = node_ref.ref();
                    way_nodes.push_back(id);
                    (*coastline_nodes)[id] = {loc.lon(), loc.lat()};
                }
            }
            if (!way_nodes.empty()) {
                (*coastline_ways)[way_nodes.front()] = way_nodes;
            }
        }
    }

   private:
    NodeMap* coastline_nodes;
    WayList* coastline_ways;
};

/**
 * @brief Class for filtering points outside of water and then creating a fmi graph from the points.
 */
class GraphCreator {
   public:
    GraphCreator() = default;
    ~GraphCreator() = default;

    /**
     * @brief Generate random points points around the globe and filter them based on the coastlines or an image.
     *
     * @param coastlines The path to the coastlines file.
     * @param num_points The number of points to generate.
     * @param output_file_prefix The prefix for the output files.
     * @param image_based_filtering Whether to use image-based filtering.
     * @param image_path The path to the image file (if using image-based filtering).
     */
    void generatePointsAndFilter(const std::string& coastlines, int num_points, const std::string& output_file_prefix,
                                 bool image_based_filtering, const std::string& image_path = "");

    /**
     * @brief Based on the filtered points, create a graph and write it to a file in the fmi format.
     *
     * @param points_file The path to the points file.
     * @param output_file_path The path to the output file.
     */
    void generateGraph(const std::string& points_file, const std::string& output_file_path);

   private:
    // The nodes and ways of the coastlines extracted from the osm.pbf file
    NodeMap m_coastline_nodes;
    WayList m_coastline_ways;

    /**
     * @brief Checks if a point is in water using image-based filtering.
     *
     * @param lat_deg Latitude of the point.
     * @param lon_deg Longitude of the point.
     * @param img_data Pointer to the image data.
     * @param img_width Width of the image.
     * @param img_height Height of the image.
     * @param img_channels Number of image channels.
     * @return True if the point is in water, false otherwise.
     */
    bool isPointInWaterImageBased(double lat_deg, double lon_deg, const unsigned char* img_data, int img_width,
                                  int img_height, int img_channels) const;

    /**
     * @brief Checks if a point is in water using polygon-based filtering.
     *
     * @param lat_deg Latitude of the point.
     * @param lon_deg Longitude of the point.
     * @param poly_vec_ways Vector of polygon ways (coastlines).
     * @param poly_boxes Vector of bounding boxes for the polygons.
     * @return True if the point is in water, false otherwise.
     */
    bool isPointInWaterPolygonBased(double lat_deg, double lon_deg,
                                    const std::vector<std::vector<labosm::Vec3>>& poly_vec_ways,
                                    const std::vector<labosm::BoudingBox>& poly_boxes) const;

    /**
     * @brief Convert latitude and longitude to pixel coordinates.
     * This uses the Equirectangular projection (Plate Carrée).
     *
     * @param lat Latitude in degrees.
     * @param lon Longitude in degrees.
     * @param img_width Width of the image in pixels.
     * @param img_height Height of the image in pixels.
     * @return A pair of pixel coordinates (x, y).
     */
    std::pair<int, int> latlon_to_pixel(double lat, double lon, int img_width, int img_height) const;

    /**
     * @brief Read points from a GeoJSON file.
     *
     * @param filename The path to the GeoJSON file.
     * @return A vector of pairs of latitude and longitude.
     */
    std::vector<std::pair<double, double>> readPointsGeoJSON(const std::string& filename);

    /**
     * @brief Write the coastlines to a GeoJSON file.
     * The coastlines are represented as LineString features.
     * @param output_file The path to the output GeoJSON file.
     * @param nodes The map of nodes (id -> coordinates).
     * @param ways The map of ways (id -> list of node ids).
     */
    void writeCoastlinesToGeojson(const std::string& output_file, const NodeMap& nodes, const WayList& ways);

    /**
     * @brief Merge the ways in the coastline data.
     * Looks at the last node of a way and checks if it is the first node of another way.
     * If so, the two ways are merged.
     */
    void merge_ways();

    /**
     * @brief Create a graph from the points.
     * Uses a KDTree to find the nearest points and create edges between them.
     * The graph is represented as a vector of vectors of edges.
     *
     * @param pointsDeg The points to create the graph from.
     * @return A vector of vectors of edges.
     */
    std::vector<std::vector<labosm::Edge>> createGraph(std::vector<std::pair<double, double>>& pointsDeg);

    /**
     * @brief Write the graph to a file in the fmi format.
     *
     * @param points The points as lat lon pairs to write to the file.
     * @param graph The graph to write to the file.
     * @param filename The path to the output file.
     */
    void writeGraphToFMI(const std::vector<std::pair<double, double>>& points,
                         const std::vector<std::vector<labosm::Edge>>& graph, const std::string& filename);

    /**
     * @brief Randomize the order of nodes using DFS.
     *
     * @param points The points to randomize.
     * @param graph The graph corresponding to the points.
     */
    void randomizeNodeOrderDFS(std::vector<std::pair<double, double>>& points,
                               std::vector<std::vector<labosm::Edge>>& graph);
};
}  // namespace labosm