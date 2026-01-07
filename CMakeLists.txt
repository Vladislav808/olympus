cmake_minimum_required(VERSION 3.16)
project(olympus LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(olympus
    src/Vector.cpp
    src/Matrix.cpp
)

target_include_directories(olympus
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

add_executable(olympus_example examples/basic.cpp)

target_link_libraries(olympus_example PRIVATE olympus)
