#include <osmium/io/any_input.hpp>
#include <osmium/handler.hpp>
#include <osmium/visitor.hpp>
#include <osmium/index/map/sparse_mem_array.hpp>
#include <osmium/handler/node_locations_for_ways.hpp>
#include <osmium/osm/location.hpp>

#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <string>

using NodeMap = std::unordered_map<uint64_t, std::pair<double, double>>;
using WayList = std::vector<std::vector<uint64_t>>;

NodeMap coastline_nodes;
WayList coastline_ways;

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
                coastline_ways.push_back(std::move(way_nodes));
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
        for (size_t i = 0; i < way.size(); ++i) {
            uint64_t id = way[i];
            auto it = nodes.find(id);
            if (it != nodes.end()) {
                const auto& [lon, lat] = it->second;
                if (i > 0) out << ",";
                out << "[" << lon << "," << lat << "]";
            }
        }
        out << R"(]},"properties":{}})";
    }

    out << "]}" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: extract_coastlines input.osm.pbf output.geojson\n";
        return 1;
    }

    const char* input_filename = argv[1];
    const char* output_filename = argv[2];

    osmium::io::Reader reader{input_filename, osmium::osm_entity_bits::node | osmium::osm_entity_bits::way};

    using index_type = osmium::index::map::SparseMemArray<osmium::unsigned_object_id_type, osmium::Location>;
    index_type index;
    osmium::handler::NodeLocationsForWays<index_type> location_handler{index};
    location_handler.ignore_errors();

    CoastlineHandler handler;

    osmium::apply(reader, location_handler, handler);
    reader.close();

    write_geojson(output_filename, coastline_nodes, coastline_ways);

    return 0;
}