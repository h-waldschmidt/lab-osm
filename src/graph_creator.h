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
using NodeMap = std::unordered_map<uint64_t, std::pair<double, double>>;
using WayList = std::unordered_map<uint64_t, std::vector<uint64_t>>;

using KDTree =
    nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3, size_t>;

class CoastlineHandler : public osmium::handler::Handler {
   public:
    CoastlineHandler(NodeMap* nodes, WayList* ways) : coastline_nodes(nodes), coastline_ways(ways) {}
    ~CoastlineHandler() = default;

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

class GraphCreator {
   public:
    GraphCreator() = default;
    ~GraphCreator() = default;

    void generatePoints(const std::string& coastlines, int num_points, const std::string& output_file_prefix,
                        bool image_based_filtering, const std::string& image_path = "");

    void generateGraph(const std::string& points_file, const std::string& output_file_path);

    std::vector<std::pair<double, double>> generatePointsOnSphere(int num_points);
    std::vector<std::pair<double, double>> readPointsGeoJSON(const std::string& filename);
    void write_geojson(const std::string& output_file, const NodeMap& nodes, const WayList& ways);
    void merge_ways();
    std::vector<std::pair<double, double>> filterOutsideWater(const std::vector<std::pair<double, double>>& points);
    std::vector<std::pair<double, double>> filterOutsideWaterImageBased(
        const std::vector<std::pair<double, double>>& points, const std::string& image_path);
    std::vector<std::vector<labosm::Edge>> createGraph(std::vector<std::pair<double, double>>& pointsDeg);
    void printGraphToFMI(const std::vector<std::pair<double, double>>& points,
                         const std::vector<std::vector<labosm::Edge>>& graph, const std::string& filename);

   private:
    NodeMap m_coastline_nodes;
    WayList m_coastline_ways;

    std::pair<int, int> latlon_to_pixel(double lat, double lon, int img_width, int img_height);
};
}  // namespace labosm