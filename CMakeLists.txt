cmake_minimum_required(VERSION 3.19)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(labosm LANGUAGES C CXX)

include(CheckIPOSupported)

check_ipo_supported(RESULT supported OUTPUT error)

if(supported AND CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not enabled: <${error}>")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto")

file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB_RECURSE HEADERS ${PROJECT_SOURCE_DIR}/src/*.h)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(OpenMP REQUIRED)
find_package(nanoflann)
find_package(Threads REQUIRED)
find_package(Protobuf REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(EXPAT REQUIRED)
find_package(Protozero REQUIRED)
find_package(Osmium REQUIRED)

add_executable(${PROJECT_NAME} ${HEADERS} ${SOURCES})

set(PROTOZERO_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third-party/protozero/include")

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/third-party/cpp-httplib)
target_include_directories(${PROJECT_NAME} PRIVATE ${OSMIUM_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${BZIP2_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${EXPAT_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${Protobuf_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROTOZERO_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX Threads::Threads ${Protobuf_LIBRARIES} ZLIB::ZLIB ${osmpbf_LIBRARIES} BZip2::BZip2 EXPAT::EXPAT nanoflann::nanoflann)
