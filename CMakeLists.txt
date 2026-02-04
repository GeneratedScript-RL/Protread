cmake_minimum_required(VERSION 3.12...4.0)
project(Protread)

find_package(OpenMP REQUIRED)

add_library(ProtreadCPP STATIC src/Protread.cpp)
set_target_properties(ProtreadCPP PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(MSVC)
    target_compile_options(ProtreadCPP PRIVATE /O2 /arch:AVX2 /openmp)
else()
    target_compile_options(ProtreadCPP PRIVATE -O3 -march=native -ffast-math -fopenmp)
endif()

find_package(Python COMPONENTS Interpreter Development REQUIRED)
add_subdirectory(pybind11)

pybind11_add_module(Protread src/binder.cpp)

target_link_libraries(Protread PRIVATE ProtreadCPP OpenMP::OpenMP_CXX)