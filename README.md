# Lab Course Algorithms on OSM Data

> [!NOTE]  
> Don't forget to clone with `--recursive`

## Building
Docker is just used for building, so we don't need to install all the dependencies.

1. `git clone --recursive https://github.com/h-waldschmidt/lab-osm.git`
2. `docker build -t labosm-builder .`
3. `docker run --rm -v ./:/app labosm-builder`

## Running
Run `./build/labosm` to see most usage examples.

Example workflow looks like this to generate graph with `6M` points and run simpleserver with just dijkstra algorithm:
1. `./build/labosm generate_points planet-coastlinespbf-cleaned.osm.pbf output/6M 6000000`
2. `./build/labosm points_to_fmi output/6M_filtered_points.geojson output/6M_graph.fmi`
3. `./build/labosm simpleserver output/6M_graph.fmi`

## Benchmarks
- Point in Polygon test with 16 Threads 6M Points: `36616160ms = 10h 10m`
  - Results in `4248964` Points in water `grep -o "Point" < filtered_points_6M.geojson | wc -l`
- On Land filtering based on image test 6M Points: `7993ms = 7,9s`
    - Results in `4157317` Points in water `grep -o "Point" < filtered_points_image_6M.geojson | wc -l`
- Generating FMI Graph with 16 Threads and `4248964` Points: `8911ms = 8,9s`
  - Results in `23117984` Edges
- Generating FMI Graph with 16 Threads and `4157317` Points: `9037ms = 9s`
  - Results in `22612494` Edges

## Code Structure

> [!NOTE]  
> Look into the header files for Documentation of most functions


