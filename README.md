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

Example workflow looks like this to generate graph with `4M` points and run simpleserver with just dijkstra algorithm:
1. `./build/labosm generate_points planet-coastlinespbf-cleaned.osm.pbf output/4M 4000000`
2. `./build/labosm points_to_fmi output/4M_filtered_points.geojson output/6M_graph.fmi`
3. `./build/labosm simpleserver output/4M_graph.fmi`

## Benchmarks
- Point in Polygon test with 16 Threads 4M Points: `35470860ms = 9h 51m`
  - Results in `4000000` Points in water `grep -o "Point" < filtered_points_4M.geojson | wc -l`
- On Land filtering based on image test 4M Points: `3155ms = 3,1s`
    - Results in `4000000` Points in water `grep -o "Point" < filtered_points_image_4M.geojson | wc -l`
- Generating FMI Graph with 16 Threads and `4000000` Points: `9874 = 9,8s`
  - Results in `21754558` Edges

### Preprocessing Time

#### Stuttgart Graph (Nodes: 1132113 Edges: 2292887)

| Heuristic       | CH  | Hub labels | Edges   |
| --------------- | --- | ---------- | ------- |
| IN_OUT          | 4s  | 84s        | 4603434 |
| EDGE_DIFFERENCE | 8s  | 57s        | 4633215 |
| WEIGHTED_COST   | 16s | 50s        | 4719853 |
| MIXED           | 11s | 39s        | 4435926 |

#### BW Graph (Nodes: 3600520 Edges: 7300290)

| Heuristic       | CH  | Hub labels | Edges    |
| --------------- | --- | ---------- | -------- |
| IN_OUT          | 18s | 631s       | 14689308 |
| EDGE_DIFFERENCE | 30s | 293s       | 14764608 |
| WEIGHTED_COST   | 60s | 266s       | 15056195 |
| MIXED           | 39s | 219s       | 14165279 |

#### Ocean Graph 1M (Nodes: 1000000 Edges: )

| Heuristic       | CH   | Hub labels | Edges |
| --------------- | ---- | ---------- | ----- |
| IN_OUT          | 665s |            |       |
| EDGE_DIFFERENCE | 471s |            |       |
| WEIGHTED_COST   | 946s |            |       |
| MIXED           | 416s |            |       |

#### Ocean Graph 4M (Nodes: 4000000 Edges: 21757192)

| Heuristic       | CH      | Hub labels | Edges     |
| --------------- | ------- | ---------- | --------- |
| IN_OUT          | 136876s |            | 110565282 |
| EDGE_DIFFERENCE | 39733s  |            | 90368552  |
| WEIGHTED_COST   | 109283s |            | 95500429  |
| MIXED           | 58858s  |            | 83502065  |


### Query Time

10000 Random queries

#### Stuttgart Graph

| Algorithm                       | Avg. Query Time (us) | Avg. Path Time (us) | Avg. PQ Pops | Avg. Label Size | Speed-up |
| ------------------------------- | -------------------- | ------------------- | ------------ | --------------- | -------- |
| Dijkstra                        | 86065.7              | 23.1336             | 585490       | N/A             | 1x       |
| CH (IN_OUT)                     | 146.305              | 133.126             | 856.228      | N/A             | 588.262x |
| Hub Labels (IN_OUT CH)          | 1.4074               | 105.107             | 0            | 92.1854         | 61152.3x |
| CH (EDGE_DIFFERENCE)            | 110.614              | 135.859             | 706.994      | N/A             | 778.075x |
| Hub Labels (EDGE_DIFFERENCE CH) | 1.3343               | 122.099             | 0            | 71.8954         | 64502.5x |
| CH (WEIGHTED_COST)              | 93.2113              | 115.63              | 634.294      | N/A             | 923.34x  |
| Hub Labels (WEIGHTED_COST CH)   | 1.4695               | 121.684             | 0            | 67.6682         | 58568x   |
| CH (MIXED)                      | 66.2676              | 108.754             | 495.195      | N/A             | 1298.76x |
| Hub Labels (MIXED CH)           | 1.3206               | 103.183             | 0            | 59.4003         | 65171.7x |

#### BW Graph

| Algorithm                       | Avg. Query Time (us) | Avg. Path Time (us) | Avg. PQ Pops | Avg. Label Size | Speed-up |
| ------------------------------- | -------------------- | ------------------- | ------------ | --------------- | -------- |
| Dijkstra                        | 314738               | 46.9453             | 1.85894e+06  | N/A             | 1x       |
| CH (IN_OUT)                     | 315.235              | 203.628             | 1733.28      | N/A             | 998.422x |
| Hub Labels (IN_OUT CH)          | 1.5703               | 181.314             | 0            | 144.997         | 200432x  |
| CH (EDGE_DIFFERENCE)            | 239.353              | 226.349             | 1328.79      | N/A             | 1314.95x |
| Hub Labels (EDGE_DIFFERENCE CH) | 1.2447               | 209.831             | 0            | 91.9075         | 252862x  |
| CH (WEIGHTED_COST)              | 216.412              | 212.72              | 1251.92      | N/A             | 1454.34x |
| Hub Labels (WEIGHTED_COST CH)   | 1.1959               | 184.915             | 0            | 87.8564         | 263181x  |
| CH (MIXED)                      | 146.293              | 200.213             | 937.949      | N/A             | 2151.43x |
| Hub Labels (MIXED CH)           | 1.225                | 210.393             | 0            | 81.281          | 256929x  |

#### Ocean Graph 1M

TODO: Rewrite the graph implementation to use coninuous array instead of array of array

#### Ocean Graph 4M
Before optimizing:

| Algorithm            | Avg. Query Time (us) | Avg. Path Time (us) | Avg. PQ Pops | Avg. Label Size | Speed-up |
| -------------------- | -------------------- | ------------------- | ------------ | --------------- | -------- |
| Dijkstra             | 1.02158e+06          | 87.8555             | 5.4125e+06   | N/A             | 1x       |
| CH (IN_OUT)          | 156558               | 403.611             | 235432       | N/A             | 6.52526x |
| CH (EDGE_DIFFERENCE) | 76031.6              | 373.679             | 159684       | N/A             | 13.4363x |
| CH (WEIGHTED_COST)   | 75230.1              | 393.461             | 146261       | N/A             | 13.5794x |
| CH (MIXED)           | 59532.7              | 364.091             | 125778       | N/A             | 17.16x   |

DFS reordering:

| Algorithm | Avg. Query Time (us) | Avg. Path Time (us) | Avg. PQ Pops | Avg. Label Size | Speed-up |
| --------- | -------------------- | ------------------- | ------------ | --------------- | -------- |
| Dijkstra  | 446868               | 35.32               | 5.00705e+06  | N/A             | 1x       |

DFS reordering and Radix PQ:
| Algorithm | Avg. Query Time (us) | Avg. Path Time (us) | Avg. PQ Pops | Avg. Label Size | Speed-up |
| --------- | -------------------- | ------------------- | ------------ | --------------- | -------- |
| Dijkstra  | 333678               | 35.3                | 5.06685e+06  | N/A             | 1x       |


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

### Dijkstra Route Examples

Here just some pictures of the routes returned by Dijkstra to kinda showcase the created edges:

![](images/Dijkstra_Route_1.png)

![](images/Dijkstra_Route_2.png)

![](images/Dijkstra_Route_3.png)

![](images/Dijkstra_Route_4.png)

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
