#include <cstdint>
#include <fstream>
#include <iostream>
#include <osmium/handler.hpp>
#include <osmium/handler/node_locations_for_ways.hpp>
#include <osmium/index/map/sparse_mem_array.hpp>
#include <osmium/io/any_input.hpp>
#include <osmium/osm/location.hpp>
#include <osmium/visitor.hpp>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using NodeMap = std::unordered_map<uint64_t, std::pair<double, double>>;
using WayList = std::unordered_map<uint64_t, std::vector<uint64_t>>;

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

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: extract_coastlines input.osm.pbf output.geojson\n";
        return 1;
    }

    const char* input_filename = argv[1];
    const char* output_filename = argv[2];

    auto start_reading = std::chrono::high_resolution_clock::now();

    osmium::io::Reader reader{input_filename, osmium::osm_entity_bits::node | osmium::osm_entity_bits::way};

    using index_type = osmium::index::map::SparseMemArray<osmium::unsigned_object_id_type, osmium::Location>;
    index_type index;
    osmium::handler::NodeLocationsForWays<index_type> location_handler{index};
    location_handler.ignore_errors();

    CoastlineHandler handler;

    osmium::apply(reader, location_handler, handler);
    reader.close();
    auto end_reading = std::chrono::high_resolution_clock::now();

    std::cout << "Reading time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_reading - start_reading).count() << " ms\n";

    auto start_merging = std::chrono::high_resolution_clock::now();
    merge_ways();
    auto end_merging = std::chrono::high_resolution_clock::now();
    std::cout << "Merging time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_merging - start_merging).count() << " ms\n";

    /*
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
    */

    write_geojson(output_filename, coastline_nodes, coastline_ways);

    return 0;
}