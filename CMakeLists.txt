cmake_minimum_required(VERSION 3.20)
project(nbody LANGUAGES CXX)

# Sources
set(SOURCES
    src/main.cpp
    src/Force/Force.cpp
    src/Integrator/Integrator.cpp
    src/Simulation/Simulation.cpp
)

add_executable(main ${SOURCES})

target_include_directories(main PRIVATE src)

# Force C++23 specifically for this executable
target_compile_features(main PRIVATE cxx_std_23)

# Prevent compiler-specific extensions (forces -std=c++23 instead of -std=gnu++23)
set_target_properties(main PROPERTIES CXX_EXTENSIONS OFF) 

# Compiler optimizations and strict warnings
if(MSVC)
    target_compile_options(main PRIVATE 
        /O2 
        /arch:AVX2
        /W4
        /WX
    )
else()
    target_compile_options(main PRIVATE 
        -O3 
        -march=native
        -std=c++23
        -Wall 
        -Wextra 
        -Wpedantic 
        -Werror
    )
endif()

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)