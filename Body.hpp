#pragma once

#include "Classes/Particle/Particle.hpp"
#include "Constants.hpp"

struct Body {
    const char* name;
    double mass;
    double x, y, z;
    double v_x, v_y, v_z;
};

// NASA Initial Conditions:
inline constexpr Body bodies[] = {
    { "Sun", 1.989e30,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0 },
    { "Mercury", 3.302e23,
        -1.478e7, -6.553e7, -3.898e6,
        36.2, -9.0, -4.0 },
    { "Venus", 4.869e24,
        -5.765e7, -9.361e7, 2.110e6,
        29.9, -18.5, -2.1 },
    { "Earth", 5.972e24,
        -2.627e7, 1.445e8, -1.049e4,
        -29.8, -5.4, 0.0 },
    { "Mars", 6.417e23,
        2.067e8, 4.500e7, -4.057e6,
        -3.9, 26.0, 0.65 },
    { "Jupiter", 1.898e27,
        5.765e8, 4.405e8, -1.493e7,
        -7.9, 10.7, 0.14 },
    { "Saturn", 5.683e26,
        1.357e9, -5.194e8, -4.480e7,
        2.9, 9.0, -0.26 },
    { "Uranus", 8.681e25,
        1.855e9, 2.233e9, -1.579e7,
        -5.2, 4.0, 0.08 },
    { "Neptune", 1.024e26,
        4.461e9, -2.705e8, -9.775e7,
        0.29, 5.5, -0.12 },
    { "Pluto", 1.303e22,
        2.595e9, -4.513e9, -2.816e8,
        4.8, 1.5, -1.6 }
};

inline void initialize_bodies( Particles &particles, std::size_t const num_bodies ) {
    for ( std::size_t i{}; i < num_bodies; ++i ) {
        particles.mass()[i] = bodies[i].mass;
        particles.pos_x()[i] = bodies[i].x * constant::KM_TO_M;
        particles.pos_y()[i] = bodies[i].y * constant::KM_TO_M;
        particles.pos_z()[i] = bodies[i].z * constant::KM_TO_M;
        particles.vel_x()[i] = bodies[i].v_x * constant::KM_TO_M;
        particles.vel_y()[i] = bodies[i].v_y * constant::KM_TO_M;
        particles.vel_z()[i] = bodies[i].v_z * constant::KM_TO_M;
        particles.acc_x()[i] = 0.0;
        particles.acc_y()[i] = 0.0;
        particles.acc_z()[i] = 0.0;
    }
}