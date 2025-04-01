# cmake-vcpkg-template
[![cpp-linter](https://github.com/h-waldschmidt/cmake-vcpkg-template/actions/workflows/lint.yml/badge.svg)](https://github.com/h-waldschmidt/cmake-vcpkg-template/actions/workflows/lint.yml)
[![codecov](https://codecov.io/gh/h-waldschmidt/cmake-vcpkg-template/graph/badge.svg?token=ECSMrHzpIn)](https://codecov.io/gh/h-waldschmidt/cmake-vcpkg-template)

Template for integrating vcpkg into a CMake project
Also includes multiple presets for building and activates compiler warnings.

Make sure `ninja` is installed for building.

## Cloning/Updating a repository which uses this template

Remember to recursively clone the repository to include `vcpkg`.

```
git clone --recursive https://github.com/h-waldschmidt/cmake-vcpkg-template.git
```

To update the `vcpkg` submodule do:

```
git submodule update --remote --merge
```

## Building

1. `cmake --preset Debug`
2. `cmake --build --preset Debug`

## Coverage

1. `cmake --preset DebugTest`
2. `cmake --build --preset DebugTest --target coverage`

## Docker

TODO:

## Presets

1. **Debug:** Compile with debug information and no optimizations.
2. **RelWithDebInfo:** Compile with debug information and optimizations.
3. **Release:** Compile without debug information and with optimizations.
4. **DebugTest** Similar to `Debug` but with an additional `${PROJECT_NAME}_test` target, which runs all defined tests. Also adds `coverage` target, which generates HTML code coverage report.
5. **ReleaseTest** Same as `DebugTest`, but with optimizations and without debug information.

## Adding Dependencies

1. Check if dependencey is available with [vcpkg](https://vcpkg.io/en/packages)
2. Add dependency to `vcpkg.json`
3. Add corresponding lines to `CMakeLists.txt` (usually `find_package` and `target_link_libraries`)
4. Rebuild
