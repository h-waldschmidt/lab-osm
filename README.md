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
- Point in Polygon test with 16 Threads: `24213776ms = 6h 43m`
- On Land filtering based on image test: `4378ms = 4s`