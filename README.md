# Lab Course Algorithms on OSM Data


## Cloning

Remember to recursively clone the repository to include `vcpkg`.

```
git clone --recursive https://github.com/h-waldschmidt/cmake-vcpkg-template.git
```

If you didn't clone recursively just run:
```
git submodule update --init
```

To update the `vcpkg` submodule do:

```
git submodule update --remote --merge
```

## Building

1. `cmake --preset Release`
2. `cmake --build --preset Release`

## Docker

TODO:

## Benchmarks
- Point in Polygon test with 16 Threads 6M Points: `36616160ms = 10h 10m`
  - Results in `4248964` Points in water `grep -o "Point" < filtered_points_6M.geojson | wc -l`
- On Land filtering based on image test 6M Points: `7993ms = 7,9s`
    - Results in `4157317` Points in water `grep -o "Point" < filtered_points_image_6M.geojson | wc -l`
- Generating FMI Graph with 16 Threads and `4248964` Points: `8911ms = 8,9s`
  - Results in `23117984` Edges
- Generating FMI Graph with 16 Threads and `4157317` Points: `9037ms = 9s`
  - Results in `22612494` Edges
