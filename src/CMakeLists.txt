find_package(quill CONFIG REQUIRED)
find_package(httplib CONFIG REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE sublibrary)

add_subdirectory(sublibrary)