cmake_minimum_required(VERSION 3.25)
project(nbody LANGUAGES CXX)

# OpenMP
find_package(OpenMP REQUIRED)

# Sources
set(SOURCES
    src/main.cpp
    src/Output/Output.cpp
    src/Force/Force.cpp
    src/Integrator/Integrator.cpp
    src/Simulation/Simulation.cpp
)

set(CORE_SOURCES
    src/Force/Force.cpp
    src/Integrator/Integrator.cpp
)

# Common compile settings
function(configure_target target_name)
    target_include_directories(${target_name} PRIVATE src)
    target_compile_features(${target_name} PRIVATE cxx_std_23)
    set_target_properties(${target_name} PROPERTIES CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(${target_name} PRIVATE /O2 /arch:AVX2 /W4 /WX)
    else()
        target_compile_options(${target_name} PRIVATE -O3 -march=native -Wall -Wextra -Wpedantic -Werror)
    endif()
    target_link_libraries(${target_name} PRIVATE OpenMP::OpenMP_CXX)
endfunction()

# Main simulation
add_executable(main ${SOURCES})
configure_target(main)

# Unit tests
add_executable(tests tests/unit_tests.cpp ${CORE_SOURCES})
configure_target(tests)

# Scaling benchmark
add_executable(benchmark tests/benchmark.cpp ${CORE_SOURCES})
configure_target(benchmark)