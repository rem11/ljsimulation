cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

set(CMAKE_CXX_FLAGS "-O3")

project(ljsimple LANGUAGES CXX)
find_package(OpenMP REQUIRED)
add_executable(ljsimple ljsimple.cpp common.cpp)
target_link_libraries(ljsimple PRIVATE OpenMP::OpenMP_CXX)
