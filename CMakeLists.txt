cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CUDA_FLAGS "-prec-div=false")

project(ljsimple LANGUAGES CXX)
find_package(OpenMP REQUIRED)
add_executable(ljsimple ljsimple.cpp common.cpp)
target_link_libraries(ljsimple PRIVATE OpenMP::OpenMP_CXX)

project(ljcuda LANGUAGES CXX CUDA)
add_executable(ljcuda ljcuda.cu common.cpp)
set_target_properties(ljcuda PROPERTIES CUDA_ARCHITECTURES "86")
