FROM debian:buster-slim AS build

RUN apt-get update; \
    apt-get install -y git zip curl unzip tar build-essential pkg-config ninja-build; \
    apt-get clean;

# install newest cmake version
ARG CMAKE_VERSION=3.28.3
RUN ARCH=$(uname -m) && \
    curl -sSL https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-${ARCH}.sh -o cmake-install.sh && \
    chmod +x cmake-install.sh && \
    ./cmake-install.sh --prefix=/usr/local --skip-license && \
    rm cmake-install.sh

COPY CMakeLists.txt /app/CMakeLists.txt
COPY *.json /app/
COPY src /app/src
COPY vcpkg /app/vcpkg

WORKDIR /app

RUN cmake --preset Release && \
    cmake --build --preset Release

FROM debian:buster-slim AS release

COPY --from=build /app/build /app/build

WORKDIR /app

ENTRYPOINT [ "./build/src/test-project" ]
