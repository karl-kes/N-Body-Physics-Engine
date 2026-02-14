cmake_minimum_required(VERSION 3.16)
project(nbody LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Sources
set(SOURCES
    src/main.cpp
    src/Force/Force.cpp
    src/Integrator/Integrator.cpp
    src/Particle/Particle.cpp
    src/Simulation/Simulation.cpp
)

add_executable(main ${SOURCES})

target_include_directories(main PRIVATE src)

# Compiler optimizations
if(MSVC)
    target_compile_options(main PRIVATE /O2 /arch:AVX2)
else()
    target_compile_options(main PRIVATE -O3 -march=native)
endif()

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)