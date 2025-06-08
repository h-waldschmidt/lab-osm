

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

#include "graph_creator.h"
#include "server.h"

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

int main(int argc, char* argv[]) {
    if (argc <= 2) {
        std::cout << "Usage ./labosm <mode> <args>\n" << "\n";
        std::cout << "Supported modes: generate_points points_to_fmi simpleserver advancedserver\n" << "\n";

        std::cout << "Example: ./labosm generate_points input.osm.pbf output_prefix num_points" << "\n";
        std::cout << "Example: ./labosm generate_points input.osm.pbf output_prefix num_points image_path\n";
        std::cout << "Generates num_points random points on the sphere, filters them based on the given coastlines and "
                     "writes them to a geojson file\n";
        std::cout << "Creates the files output_prefix_coastlines.geojson and output_prefix_filtered_points.geojson\n";
        std::cout << "If image_path is given, the points are filtered based on the image\n" << "\n";

        std::cout << "Example: ./labosm points_to_fmi filtered_points.geojson output.fmi\n" << "\n";
        std::cout << "Converts the filtered points to a fmi file\n" << "\n";

        // TODO: Take fmi file as input, generate contraction graph and write it to a fmi file
        // TODO: Support contraction graph as input for server
        // TODO: Add memory usage values and number of PQ pops to server
        // TODO: Add a benchmark function that tests all the routing algorithms with 1000 random point to point queries

        std::cout << "Example: ./labosm simpleserver input.fmi" << "\n";
        std::cout << "Simple Server only supports dijkstra for routing\n" << "\n";
        std::cout << "Example: ./labosm advancedserver input.fmi" << "\n";
        std::cerr << "Advanced Server supports CH and Hub Labeling for routing" << "\n";
        return 1;
    }

    if (argv[1] == std::string("simpleserver")) {
        if (argc != 3) {
            std::cerr << "Usage: ./labosm simpleserver input.fmi" << "\n";
            return 1;
        }
        labosm::server::simpleServer(argv[2]);
        return 0;
    } else if (argv[1] == std::string("advancedserver")) {
        if (argc != 3) {
            std::cerr << "Usage: ./labosm advancedserver input.fmi" << "\n";
            return 1;
        }
        labosm::server::advancedServer(argv[2]);
        return 0;
    } else if (argv[1] == std::string("generate_points")) {
        if (argc < 5 || argc > 6) {
            std::cerr << "Usage: ./labosm generate_points input.osm.pbf output_prefix num_points [image_path]" << "\n";
            return 1;
        }
        labosm::GraphCreator graph_creator;
        std::string image_path = (argc == 6) ? argv[5] : "";
        graph_creator.generatePointsAndFilter(argv[2], std::stoi(argv[4]), argv[3], image_path != "", image_path);
        return 0;
    } else if (argv[1] == std::string("points_to_fmi")) {
        if (argc != 4) {
            std::cerr << "Usage: ./labosm points_to_fmi filtered_points.geojson output.fmi" << "\n";
            return 1;
        }
        labosm::GraphCreator graph_creator;
        graph_creator.generateGraph(argv[2], argv[3]);
        return 0;
    } else {
        std::cerr << "Unknown mode: " << argv[1] << "\n";
        return 1;
    }

    return 0;
}
