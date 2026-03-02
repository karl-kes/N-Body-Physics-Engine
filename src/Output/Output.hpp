#pragma once

#include "Particle/Particle.hpp"

#include <fstream>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

/*
    Binary format:
    Header
    uint64_t  num_bodies
    For each body:
        char[32] name  (null-padded)

    Per frame
    uint64_t  step
    double    time_s
    For each body:
        double x, y, z, vx, vy, vz   (meters, m/s)
*/

class Binary_Output {
private:
    std::ofstream file_;
    std::size_t num_bodies_;
    std::vector<double> buffer_;

public:
    explicit Binary_Output(
        std::string const &path, 
        std::vector<std::string> const &names, 
        std::size_t const num
    );

    void write(
        Particles const &particles,
        std::size_t const num,
        std::size_t const step,
        double const time
    );
};