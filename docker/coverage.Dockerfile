FROM debian:buster-slim

RUN apt-get update; \
    apt-get install -y git zip curl unzip tar build-essential pkg-config ninja-build gcovr; \
    apt-get clean;

# install newest cmake version
ARG CMAKE_VERSION=3.28.3
RUN ARCH=$(uname -m) && \
    curl -sSL https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-${ARCH}.sh -o cmake-install.sh && \
    chmod +x cmake-install.sh && \
    ./cmake-install.sh --prefix=/usr/local --skip-license && \
    rm cmake-install.sh

VOLUME [ "/app" ]
WORKDIR /app

ENTRYPOINT cmake --preset DebugTest && \
    cmake --build --preset DebugTest && \
    ./build/src/test-project_test && \
    gcovr --xml -e /app/build/vcpkg_installed -e '.*main.cpp' -e '.*_test.cpp' > build/coverage.xml

