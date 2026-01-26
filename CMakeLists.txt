cmake_minimum_required(VERSION 3.10)
project(ThermoneutronicsSimulation)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -Wall")

set(CMAKE_BUILD_TYPE Release)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set(SOURCES src/config.cpp src/material.cpp src/mesh.cpp src/scheme.cpp src/solver.cpp src/test_modes.cpp src/main.cpp)

add_executable(simulation ${SOURCES})
target_link_libraries(simulation Eigen3::Eigen)



