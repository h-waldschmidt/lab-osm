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

For image based filtering download [this image from NASA](https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73963/gebco_08_rev_bath_21600x10800.png).

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

## Images

### 2000 Points filtered with Coastlines

Unfiltered Points:

|          Globe           |          Map           |
| :----------------------: | :--------------------: |
| ![](images/2K_Globe.png) | ![](images/2K_Map.png) |

It can be seen that the distribution is uniform and doesn't cluster at the poles.

Filtered Points:

![](images/2K_Filtered.png)

### 2000 Points filtered with Image

Unfiltered Points:

|             Globe              |             Map              |
| :----------------------------: | :--------------------------: |
| ![](images/2K_Image_Globe.png) | ![](images/2K_Image_Map.png) |

It can be seen that the distribution is uniform and doesn't cluster at the poles.

Filtered Points:

![](images/2K_Image_Filtered.png)

## Code Structure

> [!NOTE]  
> Look into the header `.h` files for Documentation of most functions

### `main.cpp`

Then entrypoint for everything (osm preprocessing, graph creation and the routing servers).

### `helper.h`

Defines a bunch of structs and small helper functions


### `graph_creator.h` and `graph_creator.cpp`

Provides all the OSM, point in water filtering and graph creation mechanisms.
Takes the `*.osm.pbf` coastlines files (and image if image based filtering is used) as input, creates `n` random points, filters out all the points that are not in water and then finally creates a graph from the filtered points and writes them to a `*.fmi` file.

### `graph.h` and `graph.cpp`

Takes an `*.fmi` files as input and creates the general graph datastructure from it.
More importantly provides functions for contraction hierarchy and hub label calculation.
Then based on these datastructer route planning using vanilla dijkstra, contraction hierarchies and hub labeling can be done.

### `server.h` and `server.cpp`

Wraps a server around the `graph.h` that is called from the frontend.
Provides two modes the `simpleserver` for just dijkstra routing and the `advancedserver` for CH and Hub-Label routing.


### `index.html` and `index.js`

Provide the frontend to display routes with leaflet.js.
