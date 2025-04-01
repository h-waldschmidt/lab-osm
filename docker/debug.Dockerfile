FROM debian:buster-slim

RUN apt-get update; \
    apt-get install -y git zip curl unzip tar build-essential pkg-config ninja-build gdbserver; \
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

# TODO: Add remote debugging capabilites, but for now just run debug build
ENTRYPOINT cmake --preset Debug && \
    cmake --build --preset Debug && \
    ./build/src/test-project
# gdbserver :8000 /app/build/src/test-project
